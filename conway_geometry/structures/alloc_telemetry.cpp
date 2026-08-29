#ifdef CONWAY_ALLOC_TELEMETRY

#include "alloc_telemetry.h"

#include <atomic>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <malloc.h>  // malloc_usable_size (dlmalloc/emscripten provide it)

// The real allocator, provided by wasm-ld's --wrap.
extern "C" {
void* __real_malloc(size_t size);
void* __real_calloc(size_t count, size_t size);
void* __real_realloc(void* ptr, size_t size);
void __real_free(void* ptr);
}

namespace {

constexpr int kSiteCount = static_cast<int>(conway::AllocSite::Count);

// ---- per-thread, active only inside an AllocTelemetryScope ----------------

struct ThreadFaceStats {
  bool active = false;
  // Which scope kind opened the active unit; every aggregate below is bucketed
  // by it so an extrusion's numbers never average into a face's.
  int kind = 0;
  uint64_t allocCalls = 0;
  uint64_t liveBytes = 0;
  uint64_t peakBytes = 0;
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

inline void onAlloc(void* ptr) {
  if (!tls.active || ptr == nullptr) {
    return;
  }
  const uint64_t size = malloc_usable_size(ptr);
  tls.allocCalls += 1;
  tls.liveBytes += size;
  if (tls.liveBytes > tls.peakBytes) {
    tls.peakBytes = tls.liveBytes;
  }
  const int site = static_cast<int>(g_currentSite);
  g_siteCounts[tls.kind][site].fetch_add(1, std::memory_order_relaxed);
  g_siteBytes[tls.kind][site].fetch_add(size, std::memory_order_relaxed);
}

inline void onFree(void* ptr) {
  if (!tls.active || ptr == nullptr) {
    return;
  }
  uint64_t size = malloc_usable_size(ptr);
  // Frees of memory allocated before the scope began can underflow the
  // in-scope live counter; clamp - we are sizing the scratch arena, and
  // pre-scope memory would not live in it.
  tls.liveBytes = (size > tls.liveBytes) ? 0 : tls.liveBytes - size;
}

// ---- process-wide aggregates, merged on scope exit -------------------------

// log2 buckets of per-scope peak bytes: bucket i covers [2^i, 2^(i+1)).
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
// Bytes still live when the scope closes - allocations that escape the unit
// (mesh growth/commits). These would NOT live in a per-unit scratch arena;
// reporting them separately decomposes scratch (arena-sizable) from commit
// (reservation-aware accumulation's job), and their ratio to peak is the
// number conway#637 asks for on each path.
std::atomic<uint64_t> g_totalEscapedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxEscapedBytes[kSiteCount] = {};

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
    "csg_operand_prep", "csg_kernel"};

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
  onFree(ptr);
  void* out = __real_realloc(ptr, size);
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
  }
}

AllocTelemetryScope::~AllocTelemetryScope() {
  if (!outermost_) {
    return;
  }
  const int kind = tls.kind;
  tls.active = false;
  g_scopes[kind].fetch_add(1, std::memory_order_relaxed);
  g_totalAllocCalls[kind].fetch_add(tls.allocCalls, std::memory_order_relaxed);
  g_totalPeakBytes[kind].fetch_add(tls.peakBytes, std::memory_order_relaxed);
  atomicMax(g_maxAllocCallsPerScope[kind], tls.allocCalls);
  atomicMax(g_maxPeakBytes[kind], tls.peakBytes);
  g_peakHistogram[kind][bucketFor(tls.peakBytes)].fetch_add(
      1, std::memory_order_relaxed);
  g_totalEscapedBytes[kind].fetch_add(tls.liveBytes, std::memory_order_relaxed);
  atomicMax(g_maxEscapedBytes[kind], tls.liveBytes);
}

namespace {

// Report one scope kind. Kinds with no scopes print nothing, so a model that
// takes only some paths gives a short report rather than a wall of zeroes --
// but note the ABSENCE of a kind here means "no scope of that kind ran", which
// on an uninstrumented call graph is indistinguishable from "no work". Only
// the presence of a nonzero line is evidence.
void DumpScopeKind(int kind) {
  const uint64_t scopes = g_scopes[kind].load();

  if (scopes == 0) {
    return;
  }

  const uint64_t allocCalls = g_totalAllocCalls[kind].load();
  const uint64_t peakTotal = g_totalPeakBytes[kind].load();
  const uint64_t escapedTotal = g_totalEscapedBytes[kind].load();

  fprintf(stderr,
          "[alloc-telemetry]   scope %-13s: units=%" PRIu64
          " allocCalls(total=%" PRIu64 " avg=%.1f max=%" PRIu64
          ") peakBytes(avg=%" PRIu64 " max=%" PRIu64 ")\n",
          kSiteNames[kind], scopes, allocCalls,
          static_cast<double>(allocCalls) / static_cast<double>(scopes),
          g_maxAllocCallsPerScope[kind].load(), peakTotal / scopes,
          g_maxPeakBytes[kind].load());
  // peak:retained is the ratio conway#635/#637 reason about. Summed over units
  // rather than per unit, so it is comparable to the whole-load figures in the
  // ledger; escapedTotal == 0 would divide by zero, so it is reported as such.
  fprintf(stderr,
          "[alloc-telemetry]     retainedBytes(avg=%" PRIu64 " max=%" PRIu64
          " total=%" PRIu64 ") peakTotal=%" PRIu64 " peak:retained=",
          escapedTotal / scopes, g_maxEscapedBytes[kind].load(), escapedTotal,
          peakTotal);

  if (escapedTotal == 0) {
    fprintf(stderr, "n/a (nothing retained)\n");
  } else {
    fprintf(stderr, "%.1fx\n",
            static_cast<double>(peakTotal) /
                static_cast<double>(escapedTotal));
  }

  for (int site = 0; site < kSiteCount; ++site) {
    const uint64_t count = g_siteCounts[kind][site].load();

    if (count == 0) {
      continue;
    }

    fprintf(stderr,
            "[alloc-telemetry]     site %-17s: %12" PRIu64
            " (%5.1f%%, %.1f/unit) %12" PRIu64 " B\n",
            kSiteNames[site], count,
            allocCalls == 0 ? 0.0
                            : 100.0 * static_cast<double>(count) /
                                  static_cast<double>(allocCalls),
            static_cast<double>(count) / static_cast<double>(scopes),
            g_siteBytes[kind][site].load());
  }

  // Cumulative histogram: what fraction of units fit an arena of 2^(i+1)?
  uint64_t running = 0;

  for (int i = 0; i < kBuckets; ++i) {
    const uint64_t count = g_peakHistogram[kind][i].load();

    if (count == 0) {
      continue;
    }

    running += count;
    fprintf(stderr,
            "[alloc-telemetry]     peak<%8" PRIu64 "KiB: %10" PRIu64
            " units (%6.2f%% cumulative)\n",
            (uint64_t(1) << (i + 1)) / 1024, count,
            100.0 * static_cast<double>(running) /
                static_cast<double>(scopes));
  }
}

}  // namespace

void DumpAllocTelemetry(const char* label) {
  uint64_t scopes = 0;

  for (int kind = 0; kind < kSiteCount; ++kind) {
    scopes += g_scopes[kind].load();
  }

  if (scopes == 0) {
    fprintf(stderr, "[alloc-telemetry] %s: no scoped units recorded\n",
            label != nullptr ? label : "");
    return;
  }

  fprintf(stderr, "[alloc-telemetry] %s: units=%" PRIu64 "\n",
          label != nullptr ? label : "", scopes);

  for (int kind = 0; kind < kSiteCount; ++kind) {
    DumpScopeKind(kind);
  }
}

AllocTagScope::AllocTagScope(AllocSite site) : previous_(g_currentSite) {
  g_currentSite = site;
}

AllocTagScope::~AllocTagScope() { g_currentSite = previous_; }

void ResetAllocTelemetry() {
  for (int kind = 0; kind < kSiteCount; ++kind) {
    g_scopes[kind].store(0);
    g_totalAllocCalls[kind].store(0);
    g_maxAllocCallsPerScope[kind].store(0);
    g_totalPeakBytes[kind].store(0);
    g_maxPeakBytes[kind].store(0);
    g_totalEscapedBytes[kind].store(0);
    g_maxEscapedBytes[kind].store(0);

    for (int i = 0; i < kBuckets; ++i) {
      g_peakHistogram[kind][i].store(0);
    }

    for (int site = 0; site < kSiteCount; ++site) {
      g_siteCounts[kind][site].store(0);
      g_siteBytes[kind][site].store(0);
    }
  }
}

}  // namespace conway

#endif  // CONWAY_ALLOC_TELEMETRY
