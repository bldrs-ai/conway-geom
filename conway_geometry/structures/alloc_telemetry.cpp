#ifdef CONWAY_ALLOC_TELEMETRY

#include "alloc_telemetry.h"

#include <atomic>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <malloc.h>  // malloc_usable_size (dlmalloc/emscripten provide it)

// Slots in the per-thread ownership table, as a power of two. The table holds
// one entry per allocation that is currently LIVE inside the open scope, not
// one per allocation the scope makes, because deletion compacts (see
// tableErase) rather than leaving tombstones -- so an allocate/free-churning
// unit occupies a constant handful of slots however many calls it makes.
//
// 2^18 entries is 6 MB per thread on a 64-bit host and 4 MB on wasm32, taken
// once from the unwrapped allocator on the first scope a thread opens. It is
// deliberately generous: when the table fills, ownership degrades to
// counted-but-unclassified for the remainder of that scope, which is reported
// but is still a worse measurement. Override at compile time if a model shows
// nonzero `unowned` in the report.
#ifndef CONWAY_ALLOC_TELEMETRY_TABLE_BITS
#define CONWAY_ALLOC_TELEMETRY_TABLE_BITS 18
#endif

// The real allocator, provided by wasm-ld's --wrap.
extern "C" {
void* __real_malloc(size_t size);
void* __real_calloc(size_t count, size_t size);
void* __real_realloc(void* ptr, size_t size);
void __real_free(void* ptr);
}

namespace {

constexpr int kSiteCount = static_cast<int>(conway::AllocSite::Count);

// ---- load-wide denominator ------------------------------------------------

// Incremented on EVERY wrapped call, inside a scope or not. This is the only
// part of the instrument that observes unscoped traffic, and it is the whole
// difference between "csg_boolean makes 43.2 M allocator calls" (a lever) and
// "csg_boolean is N % of this load's allocator calls" (a ranking). conway#637
// published a share against the in-scope subtotal and had to retract it.
std::atomic<uint64_t> g_loadAllocCalls{0};
// Calls that returned null. `g_loadAllocCalls` counts attempts, so successful
// allocations are calls - failed, and a share taken against the denominator
// uses that difference. Kept as its own term rather than folded away because
// "which paths were allocating when the heap ran out" is exactly the question
// a failure census answers, and this build can produce failures: every wasm
// target links -s ABORTING_MALLOC=0 (genie.lua).
std::atomic<uint64_t> g_loadAllocFailed{0};
std::atomic<uint64_t> g_loadAllocBytes{0};
std::atomic<uint64_t> g_loadFreeCalls{0};
std::atomic<uint64_t> g_loadFreeBytes{0};

// ---- per-thread ownership table -------------------------------------------

constexpr uint32_t kTableBits = CONWAY_ALLOC_TELEMETRY_TABLE_BITS;
constexpr uint32_t kTableSlots = uint32_t(1) << kTableBits;
constexpr uint32_t kTableMask = kTableSlots - 1;

// Linear probing degrades sharply as a table approaches full, and both probe
// loops below terminate by reaching a non-current slot. Capping live entries
// at 7/8 of capacity keeps at least kTableSlots/8 slots non-current at all
// times, which makes that termination a property of the code rather than of
// the workload, and turns "the table is full" into a definite, counted event.
constexpr uint32_t kTableLoadLimit = (kTableSlots / 8) * 7;

struct OwnedEntry {
  void* ptr;
  uint32_t bytes;
  uint32_t site;
  // Which scope activation owns this slot. A slot whose epoch is not the
  // thread's current epoch is free, which is what makes "clear the table"
  // O(1) at scope entry -- essential when the instrument opens 29,540 scopes
  // in one load. Epoch 0 is reserved for "free", so epochs start at 1.
  uint32_t epoch;
};

thread_local OwnedEntry* tlsTable = nullptr;
thread_local uint32_t tlsEpoch = 0;
thread_local uint32_t tlsTableUsed = 0;

inline uint32_t slotFor(const void* ptr) {
  // Fibonacci hashing over the pointer with the always-identical alignment
  // bits shifted out first; without the shift every key lands in one eighth
  // of the table.
  const uint64_t key =
      static_cast<uint64_t>(reinterpret_cast<uintptr_t>(ptr)) >> 4;

  return static_cast<uint32_t>((key * 0x9E3779B97F4A7C15ull) >>
                               (64 - kTableBits));
}

/** Open a fresh epoch: every slot written under an earlier one becomes free.
 *  Allocates the table on first use for this thread. */
inline void tableBeginEpoch() {
  if (tlsTable == nullptr) {
    // Through the UNWRAPPED allocator on purpose. Routing this through
    // malloc() would re-enter __wrap_malloc -> onAlloc -> tableBeginEpoch.
    // Calling __real_malloc from the wrapper is sequential, not nested: the
    // allocator call it is interposing has already returned.
    tlsTable = static_cast<OwnedEntry*>(
        __real_malloc(sizeof(OwnedEntry) * static_cast<size_t>(kTableSlots)));

    // The wasm build links -s ABORTING_MALLOC=0 (genie.lua), so a null return
    // is a state this build can actually produce. Every allocation the scope
    // then sees is counted unowned, which the report shows.
    if (tlsTable == nullptr) {
      return;
    }

    memset(tlsTable, 0, sizeof(OwnedEntry) * static_cast<size_t>(kTableSlots));
    tlsEpoch = 0;
  }

  tlsTableUsed = 0;
  tlsEpoch += 1;

  // On wrap every stale slot would read as current again, so this is the one
  // point where the table has to be cleared the expensive way.
  if (tlsEpoch == 0) {
    memset(tlsTable, 0, sizeof(OwnedEntry) * static_cast<size_t>(kTableSlots));
    tlsEpoch = 1;
  }
}

/** Record `ptr` as owned by the open scope. False when the table refused it,
 *  in which case the caller counts the allocation as unowned. */
inline bool tableInsert(void* ptr, uint32_t bytes, uint32_t site) {
  if (tlsTable == nullptr || tlsTableUsed >= kTableLoadLimit) {
    return false;
  }

  uint32_t at = slotFor(ptr);

  // A live pointer is in the table at most once: it can only be inserted
  // while a scope is open, and the free (or the realloc's free half) that
  // ends its life erases it first.
  while (tlsTable[at].epoch == tlsEpoch) {
    at = (at + 1) & kTableMask;
  }

  tlsTable[at] = OwnedEntry{ptr, bytes, site, tlsEpoch};
  tlsTableUsed += 1;

  return true;
}

/** Slot index of `ptr` if the open scope allocated it, else -1. */
inline int32_t tableFind(const void* ptr) {
  if (tlsTable == nullptr) {
    return -1;
  }

  uint32_t at = slotFor(ptr);

  while (tlsTable[at].epoch == tlsEpoch) {
    if (tlsTable[at].ptr == ptr) {
      return static_cast<int32_t>(at);
    }

    at = (at + 1) & kTableMask;
  }

  return -1;
}

/** Remove the entry at `hole`, closing the probe chain behind it.
 *
 *  Backward-shift deletion rather than tombstones, because tombstones would
 *  make occupancy grow with the number of allocations a scope MAKES instead of
 *  the number it holds live -- and csg_boolean makes ~6,700 allocations per
 *  composition while holding far fewer at once. With compaction the table is
 *  sized by live depth, which is what makes 2^18 slots a generous bound rather
 *  than a marginal one. */
inline void tableErase(uint32_t hole) {
  uint32_t scan = hole;

  for (;;) {
    scan = (scan + 1) & kTableMask;

    if (tlsTable[scan].epoch != tlsEpoch) {
      break;
    }

    const uint32_t home = slotFor(tlsTable[scan].ptr);

    // Move an entry down into the hole only when its home slot is NOT
    // cyclically inside (hole, scan]. Moving one whose home lies in that range
    // would place it before its own home, and a later lookup starting at that
    // home would walk past it.
    const bool mustStay = hole < scan ? (home > hole && home <= scan)
                                      : (home > hole || home <= scan);

    if (mustStay) {
      continue;
    }

    tlsTable[hole] = tlsTable[scan];
    hole = scan;
  }

  tlsTable[hole].ptr = nullptr;
  tlsTable[hole].epoch = 0;
  tlsTableUsed -= 1;
}

// ---- per-thread, active only inside an AllocTelemetryScope ----------------

struct ThreadFaceStats {
  bool active = false;
  // Which scope kind opened the active unit; every aggregate below is bucketed
  // by it so an extrusion's numbers never average into a face's.
  int kind = 0;
  uint64_t allocCalls = 0;
  // Bytes of in-scope allocations the ownership table accepted, minus those
  // already freed in scope. Ownership is what makes this a real quantity: only
  // a free of a pointer THIS scope activation handed out is subtracted, so it
  // cannot be driven down by a free of pre-scope memory and never needs a
  // clamp. At scope close it is exactly the escaped (still-live) bytes.
  uint64_t liveBytes = 0;
  uint64_t peakBytes = 0;
  // Every byte handed out inside the unit, owned or not, whether or not it was
  // freed again. Total in-scope allocation VOLUME -- a superset of anything an
  // arena could hold, since it includes allocations that must outlive the
  // unit. diedBytes is the arena-eligible subset.
  uint64_t cumulativeBytes = 0;
  // In-scope allocation calls that returned null. Deliberately NOT part of
  // allocCalls: that stays the census of successful allocations, which is what
  // the ownership identity (allocCalls == ownedAllocs + unownedAllocs) and
  // every byte column rest on. A failure allocates nothing, so it belongs in
  // its own column rather than inflating one that means bytes exist.
  uint64_t failedAllocs = 0;
  uint64_t ownedAllocs = 0;
  // Allocations the table refused (it was at its load limit). Counted in
  // allocCalls, cumulativeBytes and the load-wide denominator, but not
  // lifetime-classifiable, so carried as an explicit residual instead of being
  // silently attributed to one class or the other.
  uint64_t unownedAllocs = 0;
  uint64_t unownedBytes = 0;
  uint64_t freeCalls = 0;
  uint64_t freedBytes = 0;
  // Allocated AND freed inside the unit. This is the arena-eligible subset: a
  // bump arena rewound at scope exit is correct for exactly these.
  uint64_t diedCalls = 0;
  uint64_t diedBytes = 0;
  // Frees, inside the unit, of pointers the unit did not allocate -- pre-scope
  // memory released by a Cleanup() or a kernel call, caller-owned vectors
  // regrown in place. The predecessor of this counter was `clampedFrees`,
  // which could only see the subset large enough to drive the live counter
  // negative; this one is a census, because ownership is now known for every
  // free rather than inferred from a magnitude comparison.
  uint64_t foreignFrees = 0;
  uint64_t foreignBytes = 0;
};

thread_local ThreadFaceStats tls;

// Active allocation-attribution site for this thread (see AllocTagScope).
thread_local conway::AllocSite g_currentSite = conway::AllocSite::Other;

// Per-(scope kind, site) in-scope allocation counts and gross bytes,
// process-wide. Bytes as well as counts because the two answer different
// questions: counts size the allocator-call load an arena removes, bytes size
// the arena. A path can dominate one and not the other.
std::atomic<uint64_t> g_siteCounts[kSiteCount][kSiteCount] = {};
std::atomic<uint64_t> g_siteBytes[kSiteCount][kSiteCount] = {};
// The same, restricted to allocations that died inside their scope. Per-site
// lifetime is what separates a persistent cache from unit output: the global
// VertexWelder (AllocSite::VertexWeld) grows capacity that is still live at
// scope close and would otherwise be read as escaped OUTPUT, which is one of
// the five retractions conway#653 exists to answer.
std::atomic<uint64_t> g_siteDiedCounts[kSiteCount][kSiteCount] = {};
std::atomic<uint64_t> g_siteDiedBytes[kSiteCount][kSiteCount] = {};

inline void onAlloc(void* ptr) {
  // Counted before the null test AND before the scope test, deliberately. The
  // scope test is skipped because this counter's whole purpose is to see the
  // traffic the scopes do not. The null test is skipped because filtering
  // failures out here would quietly turn a census of allocation CALLS into a
  // census of successful ones, understating exactly the paths that are
  // allocating hardest when the heap is under pressure -- and null is a state
  // this build produces rather than a theoretical one (-s ABORTING_MALLOC=0
  // on every wasm target, and dlmalloc returns 0 for a request at or above
  // MAX_REQUEST).
  g_loadAllocCalls.fetch_add(1, std::memory_order_relaxed);

  if (ptr == nullptr) {
    g_loadAllocFailed.fetch_add(1, std::memory_order_relaxed);

    // Attributed to the open scope too: "which path was allocating when the
    // allocator started refusing" is the question a failure census exists to
    // answer, and a per-kind total is how it gets answered.
    if (tls.active) {
      tls.failedAllocs += 1;
    }

    return;
  }

  const uint64_t size = malloc_usable_size(ptr);

  g_loadAllocBytes.fetch_add(size, std::memory_order_relaxed);

  if (!tls.active) {
    return;
  }

  tls.allocCalls += 1;
  tls.cumulativeBytes += size;

  const int site = static_cast<int>(g_currentSite);

  g_siteCounts[tls.kind][site].fetch_add(1, std::memory_order_relaxed);
  g_siteBytes[tls.kind][site].fetch_add(size, std::memory_order_relaxed);

  // The slot records the size in 32 bits. wasm32 cannot produce a single
  // allocation above 4 GiB, but the native test harness is 64-bit, so an
  // oversized block declines ownership rather than truncating.
  if (size > 0xFFFFFFFFull ||
      !tableInsert(ptr, static_cast<uint32_t>(size),
                   static_cast<uint32_t>(site))) {
    tls.unownedAllocs += 1;
    tls.unownedBytes += size;
    return;
  }

  tls.ownedAllocs += 1;
  tls.liveBytes += size;

  if (tls.liveBytes > tls.peakBytes) {
    tls.peakBytes = tls.liveBytes;
  }
}

/** Account for the release of `ptr`, whose usable size the caller has already
 *  measured.
 *
 *  Split out for `__wrap_realloc`, which has to read the size while the block
 *  is still valid but must not commit the accounting until it knows whether
 *  the realloc succeeded. */
inline void onFreeSized(void* ptr, uint64_t size) {
  if (ptr == nullptr) {
    return;
  }

  g_loadFreeCalls.fetch_add(1, std::memory_order_relaxed);
  g_loadFreeBytes.fetch_add(size, std::memory_order_relaxed);

  if (!tls.active) {
    return;
  }

  tls.freeCalls += 1;
  tls.freedBytes += size;

  const int32_t slot = tableFind(ptr);

  if (slot < 0) {
    // Not ours. Nothing is subtracted from liveBytes -- that subtraction, done
    // unconditionally, is what corrupted every peak and retained figure this
    // instrument published before conway#653.
    tls.foreignFrees += 1;
    tls.foreignBytes += size;
    return;
  }

  const uint64_t owned = tlsTable[slot].bytes;
  const int site = static_cast<int>(tlsTable[slot].site);

  tableErase(static_cast<uint32_t>(slot));

  tls.liveBytes -= owned;
  tls.diedCalls += 1;
  tls.diedBytes += owned;
  g_siteDiedCounts[tls.kind][site].fetch_add(1, std::memory_order_relaxed);
  g_siteDiedBytes[tls.kind][site].fetch_add(owned, std::memory_order_relaxed);
}

inline void onFree(void* ptr) {
  onFreeSized(ptr, ptr != nullptr ? malloc_usable_size(ptr) : 0);
}

// ---- process-wide aggregates, merged on scope exit -------------------------

// log2 buckets of per-scope byte totals: bucket i covers [2^i, 2^(i+1)).
// 48 buckets reaches 256 TiB, comfortably beyond any wasm32 value.
constexpr int kBuckets = 48;

// All of these are per scope kind. A single blended set was fine while the
// only scope was a BREP face; it stops being fine the moment two call graphs
// with different natural units (a face, a whole extruded solid) are both
// instrumented, because "avg per face" would then average over a population
// that is not faces.
std::atomic<uint64_t> g_scopes[kSiteCount] = {};
std::atomic<uint64_t> g_totalAllocCalls[kSiteCount] = {};
std::atomic<uint64_t> g_maxAllocCallsPerScope[kSiteCount] = {};
std::atomic<uint64_t> g_totalPeakBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxPeakBytes[kSiteCount] = {};
std::atomic<uint64_t> g_peakHistogram[kSiteCount][kBuckets] = {};
// Bytes still live when the scope closes: allocations it made and did not
// free. Ownership-tracked, so this is now what the name says -- it was not
// before, when any free inside the unit reduced it.
std::atomic<uint64_t> g_totalEscapedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxEscapedBytes[kSiteCount] = {};
// Allocated and freed inside the unit: the arena-eligible subset, and its own
// histogram. This is the distribution that sizes a bump arena. The live peak
// is the wrong quantity (an arena does not reuse freed space before rewind)
// and cumulative volume is a superset (it includes what must escape); both
// were published as arena-sizing evidence and both were retracted.
std::atomic<uint64_t> g_totalDiedCalls[kSiteCount] = {};
std::atomic<uint64_t> g_totalDiedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxDiedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_diedHistogram[kSiteCount][kBuckets] = {};
// Cumulative (not live-peak) bytes per unit, and its own histogram.
std::atomic<uint64_t> g_totalCumulativeBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxCumulativeBytes[kSiteCount] = {};
std::atomic<uint64_t> g_cumulativeHistogram[kSiteCount][kBuckets] = {};
std::atomic<uint64_t> g_totalFailedAllocs[kSiteCount] = {};
std::atomic<uint64_t> g_totalOwnedAllocs[kSiteCount] = {};
std::atomic<uint64_t> g_totalOwnedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_totalUnownedAllocs[kSiteCount] = {};
std::atomic<uint64_t> g_totalUnownedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_totalFreeCalls[kSiteCount] = {};
std::atomic<uint64_t> g_totalFreedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_totalForeignFrees[kSiteCount] = {};
std::atomic<uint64_t> g_totalForeignBytes[kSiteCount] = {};

inline int bucketFor(uint64_t bytes) {
  int b = 0;
  while (bytes > 1 && b < kBuckets - 1) {
    bytes >>= 1;
    b += 1;
  }
  return b;
}

inline void atomicMax(std::atomic<uint64_t>& target, uint64_t value) {
  uint64_t cur = target.load(std::memory_order_relaxed);
  while (cur < value &&
         !target.compare_exchange_weak(cur, value, std::memory_order_relaxed)) {
  }
}

// One name per AllocSite. The static_assert below is what stops an added
// enumerator from becoming a silent out-of-bounds read in the report loops.
const char* const kSiteNames[] = {
    "other",         "earcut",         "cdt",
    "surface_eval",  "nurbs_inverse",  "tri_bounds",
    "tri_bspline",   "tri_cylinder",   "tri_sphere",
    "tri_toroidal",  "tri_conical",    "tri_revolution",
    "tri_extrusion", "advanced_face",  "extrude_solid",
    "sweep_solid",   "csg_boolean",    "extrude_cap",
    "csg_operand_prep", "csg_kernel", "vertex_weld"};

static_assert(sizeof(kSiteNames) / sizeof(kSiteNames[0]) == kSiteCount,
              "kSiteNames must have exactly one entry per AllocSite");

}  // namespace

// ---- the wrappers -----------------------------------------------------------

extern "C" {

void* __wrap_malloc(size_t size) {
  void* ptr = __real_malloc(size);
  onAlloc(ptr);
  return ptr;
}

void* __wrap_calloc(size_t count, size_t size) {
  void* ptr = __real_calloc(count, size);
  onAlloc(ptr);
  return ptr;
}

void* __wrap_realloc(void* ptr, size_t size) {
  // The old block's size must be read BEFORE the call -- afterwards the block
  // may have been freed or moved, and malloc_usable_size on it would be a
  // use-after-free read -- but the ACCOUNTING cannot be applied until the
  // outcome is known.
  const uint64_t oldSize = ptr != nullptr ? malloc_usable_size(ptr) : 0;

  void* out = __real_realloc(ptr, size);

  // Two different things return null here, and they need opposite handling.
  //
  //   - FAILURE, when a size was actually requested. realloc's contract
  //     leaves the ORIGINAL block allocated, so its ownership entry has to
  //     survive. An earlier revision ran the free accounting unconditionally
  //     and up front, which recorded a still-live block as died-in-scope,
  //     understated escaped by its size, and left its eventual in-scope free
  //     to be misread as foreign. Verified in both allocators this code runs
  //     against, rather than assumed: emscripten's dlrealloc takes the
  //     `bytes >= MAX_REQUEST` branch and returns 0 with oldmem untouched
  //     (system/lib/dlmalloc.c:5279), and glibc's realloc(p, SIZE_MAX)
  //     returns null with p still usable and its contents intact.
  //   - A SIZE-ZERO request the allocator chose to treat as a free. glibc
  //     does that and returns null, so the block really is gone and the
  //     normal path below is right. Emscripten's dlmalloc does NOT --
  //     REALLOC_ZERO_BYTES_FREES is left undefined there, so it reallocs down
  //     to a minimum chunk and returns non-null. Only the native harness
  //     reaches this branch.
  if (out == nullptr && ptr != nullptr && size != 0) {
    // Counts the attempt and its failure, and touches nothing else.
    onAlloc(out);
    return out;
  }

  onFreeSized(ptr, oldSize);
  onAlloc(out);
  return out;
}

void __wrap_free(void* ptr) {
  onFree(ptr);
  __real_free(ptr);
}

}  // extern "C"

namespace conway {

AllocTelemetryScope::AllocTelemetryScope(AllocSite kind) {
  if (!tls.active) {
    outermost_ = true;
    tls.active = true;
    tls.kind = static_cast<int>(kind);
    tls.allocCalls = 0;
    tls.liveBytes = 0;
    tls.peakBytes = 0;
    tls.cumulativeBytes = 0;
    tls.failedAllocs = 0;
    tls.ownedAllocs = 0;
    tls.unownedAllocs = 0;
    tls.unownedBytes = 0;
    tls.freeCalls = 0;
    tls.freedBytes = 0;
    tls.diedCalls = 0;
    tls.diedBytes = 0;
    tls.foreignFrees = 0;
    tls.foreignBytes = 0;
    // Must come after tls.active, and before any allocation can be attributed
    // to this scope: it is what makes every slot from the previous scope free.
    tableBeginEpoch();
  }
}

AllocTelemetryScope::~AllocTelemetryScope() {
  if (!outermost_) {
    return;
  }

  const int kind = tls.kind;

  // Cleared first so nothing the merge itself allocates is attributed to the
  // unit that just closed.
  tls.active = false;

  // liveBytes is, by the ownership invariant, exactly the bytes this scope
  // allocated and did not free: ownedBytes == diedBytes + escapedBytes.
  const uint64_t escaped = tls.liveBytes;

  g_scopes[kind].fetch_add(1, std::memory_order_relaxed);
  g_totalAllocCalls[kind].fetch_add(tls.allocCalls, std::memory_order_relaxed);
  g_totalPeakBytes[kind].fetch_add(tls.peakBytes, std::memory_order_relaxed);
  atomicMax(g_maxAllocCallsPerScope[kind], tls.allocCalls);
  atomicMax(g_maxPeakBytes[kind], tls.peakBytes);
  g_peakHistogram[kind][bucketFor(tls.peakBytes)].fetch_add(
      1, std::memory_order_relaxed);
  g_totalEscapedBytes[kind].fetch_add(escaped, std::memory_order_relaxed);
  atomicMax(g_maxEscapedBytes[kind], escaped);
  g_totalDiedCalls[kind].fetch_add(tls.diedCalls, std::memory_order_relaxed);
  g_totalDiedBytes[kind].fetch_add(tls.diedBytes, std::memory_order_relaxed);
  atomicMax(g_maxDiedBytes[kind], tls.diedBytes);
  g_diedHistogram[kind][bucketFor(tls.diedBytes)].fetch_add(
      1, std::memory_order_relaxed);
  g_totalCumulativeBytes[kind].fetch_add(
      tls.cumulativeBytes, std::memory_order_relaxed);
  atomicMax(g_maxCumulativeBytes[kind], tls.cumulativeBytes);
  g_cumulativeHistogram[kind][bucketFor(tls.cumulativeBytes)].fetch_add(
      1, std::memory_order_relaxed);
  g_totalFailedAllocs[kind].fetch_add(tls.failedAllocs,
                                      std::memory_order_relaxed);
  g_totalOwnedAllocs[kind].fetch_add(tls.ownedAllocs,
                                     std::memory_order_relaxed);
  g_totalOwnedBytes[kind].fetch_add(tls.diedBytes + escaped,
                                    std::memory_order_relaxed);
  g_totalUnownedAllocs[kind].fetch_add(tls.unownedAllocs,
                                       std::memory_order_relaxed);
  g_totalUnownedBytes[kind].fetch_add(tls.unownedBytes,
                                      std::memory_order_relaxed);
  g_totalFreeCalls[kind].fetch_add(tls.freeCalls, std::memory_order_relaxed);
  g_totalFreedBytes[kind].fetch_add(tls.freedBytes, std::memory_order_relaxed);
  g_totalForeignFrees[kind].fetch_add(
      tls.foreignFrees, std::memory_order_relaxed);
  g_totalForeignBytes[kind].fetch_add(
      tls.foreignBytes, std::memory_order_relaxed);
}

namespace {

/** Print one log2 histogram. `what` labels the quantity, `note` says what the
 *  distribution can be used for -- the two byte distributions here answer
 *  different questions and one of them has been misread twice. */
void DumpHistogram(const std::atomic<uint64_t>* histogram, uint64_t scopes,
                   const char* what, const char* note) {
  uint64_t running = 0;

  for (int i = 0; i < kBuckets; ++i) {
    const uint64_t count = histogram[i].load();

    if (count == 0) {
      continue;
    }

    running += count;
    fprintf(stderr,
            "[alloc-telemetry]     %s<%8" PRIu64 "KiB: %10" PRIu64
            " units (%6.2f%% cumulative)%s\n",
            what, (uint64_t(1) << (i + 1)) / 1024, count,
            100.0 * static_cast<double>(running) /
                static_cast<double>(scopes),
            note);
  }
}

// Report one scope kind. Kinds with no scopes print nothing, so a model that
// takes only some paths gives a short report rather than a wall of zeroes --
// but note the ABSENCE of a kind here means "no scope of that kind ran", which
// on an uninstrumented call graph is indistinguishable from "no work". Only
// the presence of a nonzero line is evidence.
void DumpScopeKind(int kind, uint64_t loadSuccessfulAllocs) {
  const uint64_t scopes = g_scopes[kind].load();

  if (scopes == 0) {
    return;
  }

  const uint64_t allocCalls = g_totalAllocCalls[kind].load();
  const uint64_t peakTotal = g_totalPeakBytes[kind].load();
  const uint64_t escapedTotal = g_totalEscapedBytes[kind].load();
  const uint64_t diedTotal = g_totalDiedBytes[kind].load();
  const uint64_t cumulativeTotal = g_totalCumulativeBytes[kind].load();
  const uint64_t unownedBytes = g_totalUnownedBytes[kind].load();
  const uint64_t unownedAllocs = g_totalUnownedAllocs[kind].load();

  // The share is successes over successes. allocCalls counts the allocations
  // this kind actually got; the denominator it is divided by is the load's
  // successful allocations, not its calls, so the two ends of the ratio mean
  // the same thing. Failed attempts are reported beside it, never folded in.
  fprintf(stderr,
          "[alloc-telemetry]   scope %-13s: units=%" PRIu64
          " allocCalls(total=%" PRIu64 " avg=%.1f max=%" PRIu64
          " failed=%" PRIu64 ") = %.3f%% of load-wide successful allocs\n",
          kSiteNames[kind], scopes, allocCalls,
          static_cast<double>(allocCalls) / static_cast<double>(scopes),
          g_maxAllocCallsPerScope[kind].load(),
          g_totalFailedAllocs[kind].load(),
          loadSuccessfulAllocs == 0
              ? 0.0
              : 100.0 * static_cast<double>(allocCalls) /
                    static_cast<double>(loadSuccessfulAllocs));

  // The lifetime split, and the line the arena question is answered from.
  // diedBytes is the arena-eligible subset: allocated inside the unit and
  // released inside it, so a bump arena rewound at scope exit would have held
  // exactly these and no more. escapedBytes is what must outlive the unit --
  // the returned Geometry, and any persistent cache that grew inside it (see
  // the per-site breakdown, which separates the two).
  fprintf(stderr,
          "[alloc-telemetry]     lifetime: died-in-scope(calls=%" PRIu64
          " bytes=%" PRIu64 " avg=%" PRIu64 " max=%" PRIu64
          ")  escaped(bytes=%" PRIu64 " avg=%" PRIu64 " max=%" PRIu64 ")\n",
          g_totalDiedCalls[kind].load(), diedTotal, diedTotal / scopes,
          g_maxDiedBytes[kind].load(), escapedTotal, escapedTotal / scopes,
          g_maxEscapedBytes[kind].load());

  // The identity that makes the split checkable by eye, and the residual that
  // is not classifiable when the ownership table filled.
  fprintf(stderr,
          "[alloc-telemetry]     volume: cumulativeBytes(avg=%" PRIu64
          " max=%" PRIu64 " total=%" PRIu64 ") = died %" PRIu64
          " + escaped %" PRIu64 " + unowned %" PRIu64 " (%" PRIu64
          " allocs)%s\n",
          cumulativeTotal / scopes, g_maxCumulativeBytes[kind].load(),
          cumulativeTotal, diedTotal, escapedTotal, unownedBytes,
          unownedAllocs,
          unownedAllocs == 0
              ? ""
              : "  <-- ownership table filled; raise"
                " CONWAY_ALLOC_TELEMETRY_TABLE_BITS and re-run");

  // Ownership-tracked, so this is the most bytes the unit's OWN live
  // allocations reached at once. It is not arena capacity -- an arena does not
  // reuse freed space before it rewinds, which is why the died line above is
  // the sizing quantity and this one is not.
  fprintf(stderr,
          "[alloc-telemetry]     livePeak(avg=%" PRIu64 " max=%" PRIu64
          " total=%" PRIu64 ")  [ownership-tracked; not arena capacity]\n",
          peakTotal / scopes, g_maxPeakBytes[kind].load(), peakTotal);

  // Frees, inside the unit, of memory the unit did not allocate. Under the old
  // instrument these were subtracted from the live counter and only the ones
  // big enough to drive it negative left a trace. Nothing is subtracted now,
  // so this is a diagnostic about the call graph (how much pre-scope memory it
  // releases), not a caveat on the columns above.
  fprintf(stderr,
          "[alloc-telemetry]     frees(calls=%" PRIu64 " bytes=%" PRIu64
          ") of which foreign(frees=%" PRIu64 " bytes=%" PRIu64
          ")  [foreign frees are counted, never subtracted]\n",
          g_totalFreeCalls[kind].load(), g_totalFreedBytes[kind].load(),
          g_totalForeignFrees[kind].load(),
          g_totalForeignBytes[kind].load());

  for (int site = 0; site < kSiteCount; ++site) {
    const uint64_t count = g_siteCounts[kind][site].load();

    if (count == 0) {
      continue;
    }

    const uint64_t bytes = g_siteBytes[kind][site].load();
    const uint64_t siteDied = g_siteDiedBytes[kind][site].load();

    fprintf(stderr,
            "[alloc-telemetry]     site %-17s: %12" PRIu64
            " (%5.1f%%, %.1f/unit) %12" PRIu64 " B  died=%" PRIu64
            " B  escaped<=%" PRIu64 " B\n",
            kSiteNames[site], count,
            allocCalls == 0 ? 0.0
                            : 100.0 * static_cast<double>(count) /
                                  static_cast<double>(allocCalls),
            static_cast<double>(count) / static_cast<double>(scopes), bytes,
            siteDied, bytes - siteDied);
  }

  // Three distributions, each answering a different question. `died<` is the
  // one that sizes a bump arena. `peak<` is the most a unit held live at once,
  // and `alloc<` is everything it allocated; neither sizes an arena, and each
  // was published as if it did once.
  DumpHistogram(g_diedHistogram[kind], scopes, "died ",
                " -- arena-eligible: allocated and freed inside the unit");
  DumpHistogram(g_peakHistogram[kind], scopes, "peak ",
                " -- live peak, not arena capacity");
  DumpHistogram(g_cumulativeHistogram[kind], scopes, "alloc",
                " -- total volume, superset of arena-eligible");
}

}  // namespace

void DumpAllocTelemetry(const char* label) {
  uint64_t scopes = 0;

  for (int kind = 0; kind < kSiteCount; ++kind) {
    scopes += g_scopes[kind].load();
  }

  const uint64_t loadAllocCalls = g_loadAllocCalls.load();
  const uint64_t loadAllocFailed = g_loadAllocFailed.load();
  const uint64_t loadSuccessfulAllocs = loadAllocCalls - loadAllocFailed;

  // The denominator, printed first because every in-scope figure below is a
  // share of it. Unlike the per-scope counters this one sees allocations made
  // on call graphs no scope covers, which is what lets an unscoped path be
  // ruled larger or smaller than an instrumented one. `calls` counts every
  // attempt and `failed` the null returns, so successes are the difference --
  // spelled out rather than left to the reader, since which of the two a
  // percentage was taken against is exactly the kind of thing this instrument
  // has had to retract before.
  //
  // "In scope or not" includes the runtime's own traffic -- module bring-up,
  // the filesystem shim, the pthread pool. That is correct (it IS allocator
  // traffic) but it is not deterministic, so on a model whose geometry work is
  // small the denominator moves between runs and any share taken against it
  // moves with it. Measured on this tree: a 25 KB IFC fixture reported
  // 2,710 / 4,260 / 3,108 calls over three runs of one binary -- shares of
  // 41.1 / 26.1 / 35.8 % for an in-scope count that was 1,113 every time --
  // while a fixture whose own work dominates (224,779 of 316,607) was
  // bit-identical across runs. Quote a share only when the load's own traffic
  // dominates, and say which run it came from.
  fprintf(stderr,
          "[alloc-telemetry] %s: load-wide allocator traffic:"
          " allocs(calls=%" PRIu64 " failed=%" PRIu64 " ok=%" PRIu64
          " bytes=%" PRIu64 ") frees(calls=%" PRIu64 " bytes=%" PRIu64
          ")  [every wrapped call, in scope or not -- includes"
          " nondeterministic runtime traffic; see the note in the source]\n",
          label != nullptr ? label : "", loadAllocCalls, loadAllocFailed,
          loadSuccessfulAllocs, g_loadAllocBytes.load(),
          g_loadFreeCalls.load(), g_loadFreeBytes.load());

  if (scopes == 0) {
    fprintf(stderr, "[alloc-telemetry] %s: no scoped units recorded\n",
            label != nullptr ? label : "");
    return;
  }

  fprintf(stderr, "[alloc-telemetry] %s: units=%" PRIu64 "\n",
          label != nullptr ? label : "", scopes);

  for (int kind = 0; kind < kSiteCount; ++kind) {
    DumpScopeKind(kind, loadSuccessfulAllocs);
  }
}

AllocTelemetryKindTotals GetAllocTelemetryKindTotals(AllocSite kind) {
  const int at = static_cast<int>(kind);

  AllocTelemetryKindTotals totals{};

  totals.scopes = g_scopes[at].load();
  totals.allocCalls = g_totalAllocCalls[at].load();
  totals.cumulativeBytes = g_totalCumulativeBytes[at].load();
  totals.failedAllocs = g_totalFailedAllocs[at].load();
  totals.ownedAllocs = g_totalOwnedAllocs[at].load();
  totals.ownedBytes = g_totalOwnedBytes[at].load();
  totals.unownedAllocs = g_totalUnownedAllocs[at].load();
  totals.unownedBytes = g_totalUnownedBytes[at].load();
  totals.diedCalls = g_totalDiedCalls[at].load();
  totals.diedBytes = g_totalDiedBytes[at].load();
  totals.maxDiedBytes = g_maxDiedBytes[at].load();
  totals.escapedBytes = g_totalEscapedBytes[at].load();
  totals.maxEscapedBytes = g_maxEscapedBytes[at].load();
  totals.peakBytesTotal = g_totalPeakBytes[at].load();
  totals.maxPeakBytes = g_maxPeakBytes[at].load();
  totals.freeCalls = g_totalFreeCalls[at].load();
  totals.freedBytes = g_totalFreedBytes[at].load();
  totals.foreignFrees = g_totalForeignFrees[at].load();
  totals.foreignBytes = g_totalForeignBytes[at].load();

  return totals;
}

AllocTelemetrySiteTotals GetAllocTelemetrySiteTotals(AllocSite kind,
                                                     AllocSite site) {
  const int atKind = static_cast<int>(kind);
  const int atSite = static_cast<int>(site);

  AllocTelemetrySiteTotals totals{};

  totals.allocs = g_siteCounts[atKind][atSite].load();
  totals.bytes = g_siteBytes[atKind][atSite].load();
  totals.diedCalls = g_siteDiedCounts[atKind][atSite].load();
  totals.diedBytes = g_siteDiedBytes[atKind][atSite].load();

  return totals;
}

AllocTelemetryLoadTotals GetAllocTelemetryLoadTotals() {
  AllocTelemetryLoadTotals totals{};

  totals.allocCalls = g_loadAllocCalls.load();
  totals.allocFailed = g_loadAllocFailed.load();
  totals.allocBytes = g_loadAllocBytes.load();
  totals.freeCalls = g_loadFreeCalls.load();
  totals.freeBytes = g_loadFreeBytes.load();

  return totals;
}

AllocTagScope::AllocTagScope(AllocSite site) : previous_(g_currentSite) {
  g_currentSite = site;
}

AllocTagScope::~AllocTagScope() { g_currentSite = previous_; }

void ResetAllocTelemetry() {
  g_loadAllocCalls.store(0);
  g_loadAllocFailed.store(0);
  g_loadAllocBytes.store(0);
  g_loadFreeCalls.store(0);
  g_loadFreeBytes.store(0);

  for (int kind = 0; kind < kSiteCount; ++kind) {
    g_scopes[kind].store(0);
    g_totalAllocCalls[kind].store(0);
    g_maxAllocCallsPerScope[kind].store(0);
    g_totalPeakBytes[kind].store(0);
    g_maxPeakBytes[kind].store(0);
    g_totalEscapedBytes[kind].store(0);
    g_maxEscapedBytes[kind].store(0);
    g_totalDiedCalls[kind].store(0);
    g_totalDiedBytes[kind].store(0);
    g_maxDiedBytes[kind].store(0);
    g_totalCumulativeBytes[kind].store(0);
    g_maxCumulativeBytes[kind].store(0);
    g_totalFailedAllocs[kind].store(0);
    g_totalOwnedAllocs[kind].store(0);
    g_totalOwnedBytes[kind].store(0);
    g_totalUnownedAllocs[kind].store(0);
    g_totalUnownedBytes[kind].store(0);
    g_totalFreeCalls[kind].store(0);
    g_totalFreedBytes[kind].store(0);
    g_totalForeignFrees[kind].store(0);
    g_totalForeignBytes[kind].store(0);

    for (int i = 0; i < kBuckets; ++i) {
      g_peakHistogram[kind][i].store(0);
      g_diedHistogram[kind][i].store(0);
      g_cumulativeHistogram[kind][i].store(0);
    }

    for (int site = 0; site < kSiteCount; ++site) {
      g_siteCounts[kind][site].store(0);
      g_siteBytes[kind][site].store(0);
      g_siteDiedCounts[kind][site].store(0);
      g_siteDiedBytes[kind][site].store(0);
    }
  }
}

}  // namespace conway

#endif  // CONWAY_ALLOC_TELEMETRY
