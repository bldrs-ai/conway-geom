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
  // Every byte handed out inside the unit, whether or not it was freed again
  // before the unit closed. This -- not peakBytes -- is what a bump arena has
  // to hold, because a bump arena's deallocation is a no-op until the scope
  // rewinds. See the arena-sizing note in DumpScopeKind.
  uint64_t cumulativeBytes = 0;
  uint64_t freeCalls = 0;
  uint64_t freedBytes = 0;
  // Frees whose size exceeded the in-scope live counter, i.e. where the clamp
  // in onFree fired. Direct evidence that pre-scope memory was freed inside
  // the unit, which is the condition under which liveBytes -- and so peak and
  // retained -- are wrong. Counted so the report can say so instead of the
  // reader having to guess.
  uint64_t clampedFrees = 0;
  uint64_t clampedBytes = 0;
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
  tls.cumulativeBytes += size;
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
  tls.freeCalls += 1;
  tls.freedBytes += size;
  // KNOWN DEFECT, measured rather than hidden. This instrument does not track
  // which pointers it handed out, so a free of memory allocated BEFORE the
  // scope began is indistinguishable from a free of in-scope memory and is
  // subtracted from the in-scope live counter all the same. That corrupts
  // liveBytes, and with it peakBytes and the retained figure at scope exit.
  // The clamp below only stops the counter going negative; it does not stop
  // the corruption, and where it fires it is erasing bytes that really were
  // live in scope.
  //
  // It bites hardest on the CSG path, where Cleanup() and the kernel free
  // operand buffers allocated before the composition began. The two clamp
  // counters exist so the report can quantify the exposure per scope kind
  // instead of leaving every byte column silently suspect. The real fix is
  // pointer-ownership tracking (only subtract a free whose pointer this scope
  // allocated); until that exists, treat the peak and retained columns as
  // unreliable on any kind reporting nonzero clamped frees.
  if (size > tls.liveBytes) {
    tls.clampedFrees += 1;
    tls.clampedBytes += size - tls.liveBytes;
    tls.liveBytes = 0;
  } else {
    tls.liveBytes -= size;
  }
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
// Cumulative (not live-peak) bytes per unit, and its own histogram. A bump
// arena never reuses freed space until rewind, so this is the distribution
// that sizes one, and it can be orders of magnitude above the live peak.
std::atomic<uint64_t> g_totalCumulativeBytes[kSiteCount] = {};
std::atomic<uint64_t> g_maxCumulativeBytes[kSiteCount] = {};
std::atomic<uint64_t> g_cumulativeHistogram[kSiteCount][kBuckets] = {};
// Exposure of the byte columns to the onFree ownership defect.
std::atomic<uint64_t> g_totalFreeCalls[kSiteCount] = {};
std::atomic<uint64_t> g_totalFreedBytes[kSiteCount] = {};
std::atomic<uint64_t> g_totalClampedFrees[kSiteCount] = {};
std::atomic<uint64_t> g_totalClampedBytes[kSiteCount] = {};

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
    tls.cumulativeBytes = 0;
    tls.freeCalls = 0;
    tls.freedBytes = 0;
    tls.clampedFrees = 0;
    tls.clampedBytes = 0;
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
  g_totalCumulativeBytes[kind].fetch_add(
      tls.cumulativeBytes, std::memory_order_relaxed);
  atomicMax(g_maxCumulativeBytes[kind], tls.cumulativeBytes);
  g_cumulativeHistogram[kind][bucketFor(tls.cumulativeBytes)].fetch_add(
      1, std::memory_order_relaxed);
  g_totalFreeCalls[kind].fetch_add(tls.freeCalls, std::memory_order_relaxed);
  g_totalFreedBytes[kind].fetch_add(tls.freedBytes, std::memory_order_relaxed);
  g_totalClampedFrees[kind].fetch_add(
      tls.clampedFrees, std::memory_order_relaxed);
  g_totalClampedBytes[kind].fetch_add(
      tls.clampedBytes, std::memory_order_relaxed);
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
  // Deliberately printed WITHOUT a peak:retained ratio. Both inputs are
  // derived from tls.liveBytes, which onFree corrupts whenever memory
  // allocated before the scope is freed inside it -- and that happens with no
  // clamp and no other signal whenever the freed block is smaller than what is
  // currently live. An earlier revision printed the ratio and graded it
  // per scope kind by clamp count; that grading was wrong, because zero clamps
  // is not an ownership all-clear. Until pointer ownership is tracked these
  // are raw diagnostics, not results, and the report says so on every kind
  // rather than inviting the reader to divide them.
  fprintf(stderr,
          "[alloc-telemetry]     peakBytesTotal=%" PRIu64
          " retainedBytes(avg=%" PRIu64 " max=%" PRIu64 " total=%" PRIu64
          ")  [NOT ownership-tracked, NOT lifetime-classified --"
          " raw diagnostics, they bound nothing]\n",
          peakTotal, escapedTotal / scopes, g_maxEscapedBytes[kind].load(),
          escapedTotal);

  // Total in-scope allocation volume per unit. NOT arena capacity, and an
  // earlier revision labelled it as such: sizing a bump arena needs the
  // arena-ELIGIBLE subset, i.e. allocations that die inside the scope, and this
  // total also contains the returned Geometry (which must outlive the unit by
  // definition) and persistent caches such as the global VertexWelder. The
  // instrument classifies neither ownership nor lifetime, so it cannot compute
  // that subset, and this figure bounds nothing.
  //
  // It is printed anyway because it is a raw, uncorrupted quantity -- computed
  // in onAlloc alone, never touching liveBytes -- and deleting it would only
  // mean the next reader recomputes it from the site table. What is removed is
  // the claim it was carrying.
  fprintf(stderr,
          "[alloc-telemetry]     cumulativeBytes(avg=%" PRIu64 " max=%" PRIu64
          " total=%" PRIu64 ")  [total in-scope volume, NOT arena-eligible --"
          " no lifetime classification]\n",
          g_totalCumulativeBytes[kind].load() / scopes,
          g_maxCumulativeBytes[kind].load(),
          g_totalCumulativeBytes[kind].load());

  // Exposure of the peak/retained figures to the ownership defect documented in
  // onFree. Nonzero clamped frees are PROOF the defect fired here. Zero is not
  // proof of the opposite: a pre-scope free smaller than what is currently live
  // is subtracted silently, no clamp, no signal. Concrete cases exist on paths
  // that report zero -- Extrude() copies its IfcProfile argument before the
  // scope opens and profile.curve.Add2d() can then reallocate that copy inside
  // it; a face scope grows caller-owned Geometry vectors the same way. So this
  // line quantifies a lower bound on the damage, and no scope kind's peak or
  // retained figure is trustworthy until ownership is tracked.
  const uint64_t clampedFrees = g_totalClampedFrees[kind].load();

  fprintf(stderr,
          "[alloc-telemetry]     frees(calls=%" PRIu64 " bytes=%" PRIu64
          ") clamped(frees=%" PRIu64 " bytes=%" PRIu64 ")%s\n",
          g_totalFreeCalls[kind].load(), g_totalFreedBytes[kind].load(),
          clampedFrees, g_totalClampedBytes[kind].load(),
          clampedFrees == 0
            ? "  (zero clamps is NOT an ownership all-clear)"
            : "  <-- ownership defect observed firing");

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

  // Two distributions. `peak<` is the most bytes live at once in a unit;
  // `alloc<` is every byte it allocated. NEITHER sizes a bump arena: the first
  // is the wrong quantity (an arena does not reuse freed space before rewind),
  // and the second is a superset that includes allocations which must outlive
  // the scope. Arena sizing needs the arena-eligible subset, which requires
  // lifetime classification this instrument does not have.
  uint64_t running = 0;

  for (int i = 0; i < kBuckets; ++i) {
    const uint64_t count = g_peakHistogram[kind][i].load();

    if (count == 0) {
      continue;
    }

    running += count;
    fprintf(stderr,
            "[alloc-telemetry]     peak <%8" PRIu64 "KiB: %10" PRIu64
            " units (%6.2f%% cumulative)\n",
            (uint64_t(1) << (i + 1)) / 1024, count,
            100.0 * static_cast<double>(running) /
                static_cast<double>(scopes));
  }

  running = 0;

  for (int i = 0; i < kBuckets; ++i) {
    const uint64_t count = g_cumulativeHistogram[kind][i].load();

    if (count == 0) {
      continue;
    }

    running += count;
    fprintf(stderr,
            "[alloc-telemetry]     alloc<%8" PRIu64 "KiB: %10" PRIu64
            " units (%6.2f%% cumulative) -- volume, not arena capacity\n",
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

  // Every figure below counts ONLY allocations made inside an
  // AllocTelemetryScope -- onAlloc returns immediately when !tls.active. There
  // is no load-wide denominator here, so these are absolute counts for the
  // instrumented paths and shares OF EACH OTHER, never shares of the load.
  // Reading them as "path X is N% of allocator traffic" is the same mistake,
  // one level up, as reading an unscoped path's silence as a measurement.
  fprintf(stderr,
          "[alloc-telemetry] %s: units=%" PRIu64
          " (in-scope allocations only; no load-wide denominator)\n",
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
    g_totalCumulativeBytes[kind].store(0);
    g_maxCumulativeBytes[kind].store(0);
    g_totalFreeCalls[kind].store(0);
    g_totalFreedBytes[kind].store(0);
    g_totalClampedFrees[kind].store(0);
    g_totalClampedBytes[kind].store(0);

    for (int i = 0; i < kBuckets; ++i) {
      g_peakHistogram[kind][i].store(0);
      g_cumulativeHistogram[kind][i].store(0);
    }

    for (int site = 0; site < kSiteCount; ++site) {
      g_siteCounts[kind][site].store(0);
      g_siteBytes[kind][site].store(0);
    }
  }
}

}  // namespace conway

#endif  // CONWAY_ALLOC_TELEMETRY
