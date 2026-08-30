#pragma once

/*
 * Opt-in allocation telemetry for the AFTP (Allocation-Free Tessellation
 * Pipeline) sizing pass — see design/new/emsdk-upgrade-scalable-allocator.md
 * (Phase 4) in bldrs-ai/conway.
 *
 * Compiled in only when CONWAY_ALLOC_TELEMETRY is defined, and the module is
 * linked with:
 *
 *   -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=free
 *
 * (genie.lua adds both when the CONWAY_ALLOC_TELEMETRY env var is set at
 * project generation time). Release builds are untouched.
 *
 * What it measures, for every instrumented unit of geometry work (an
 * AllocTelemetryScope): allocator calls, and -- since conway#653 -- three
 * things it could not measure before, without which no byte figure it printed
 * bounded anything.
 *
 *   1. OWNERSHIP. Every allocation made inside a scope is recorded in a
 *      per-thread pointer table, so a free is subtracted from the in-scope
 *      live counter only when this scope activation allocated that pointer.
 *      Before this, `onFree` subtracted every free that happened inside a
 *      scope, including frees of memory allocated before it began, and clamped
 *      at zero -- corrupting the live counter and with it the peak and the
 *      retained figure. Frees the scope does not own are now counted
 *      separately (`foreignFrees`) instead of being subtracted.
 *   2. LIFETIME. At scope close, the bytes it allocated split exactly into
 *      DIED-IN-SCOPE (allocated and freed inside the unit) and ESCAPED (still
 *      live when the unit closed). Only the first is arena-eligible: a bump
 *      arena rewound at scope exit is correct only for allocations that die
 *      inside the scope, so `diedBytes` -- not the live peak, and not total
 *      in-scope volume -- is the distribution that sizes one.
 *   3. A LOAD-WIDE DENOMINATOR. Counters incremented on EVERY wrapped
 *      allocation and free, in scope or not, and before the null test as well
 *      as before the scope test -- so it is a census of allocation CALLS, with
 *      the failures carried as their own term rather than filtered out.
 *      In-scope counts alone identify a lever; only a denominator ranks one.
 *
 * Counters are thread-local during collection and merged into process-wide
 * atomics on scope exit, so the instrument works identically on the MT build.
 * The ownership table is thread-local too: a pointer allocated on one thread
 * and freed inside a scope on another is foreign to that scope, which is the
 * correct answer -- that scope did not allocate it.
 *
 * Two independent axes, and confusing them is the trap this instrument fell
 * into once already (conway#637):
 *
 *   - The SCOPE KIND is which call graph opened the unit -- an advanced-BREP
 *     face, one solid extrusion, one swept solid, one boolean composition.
 *     Every aggregate (scopes, alloc calls, peak, died, escaped, histograms)
 *     is bucketed by it, so paths with different natural units do not blend
 *     into one meaningless average.
 *   - The SITE is which sub-step inside the unit made a given allocation
 *     (AllocTagScope). Sites are attributed within their enclosing scope kind,
 *     and each site's allocations are lifetime-classified too -- which is what
 *     lets a persistent cache such as the global VertexWelder be subtracted
 *     from a kind's escaped bytes instead of being read as unit output.
 *
 * The per-scope figures are still only ever collected inside an
 * AllocTelemetryScope: those wrappers are inert otherwise (the load-wide
 * counters in 3 are not). So a call graph with no scope on it reports no
 * SCOPE, which reads identically to a call graph that allocates nothing. Until
 * conway#639 audited it, Extrude() and the CSG/boolean path had no scope, and
 * the resulting "zero scoped faces on extrusion models" was read as a
 * measurement of those paths when it was only a statement about where the
 * instrument was. Before believing any null from this instrument, check that a
 * scope actually wraps the path in question.
 */

#include <cstddef>
#include <cstdint>

namespace conway {

#ifdef CONWAY_ALLOC_TELEMETRY

/** Coarse callsite buckets for attributing where in-scope allocations happen.
 *  Extend as needed; keep Count last. */
enum class AllocSite {
  Other = 0,
  Earcut,
  Cdt,
  SurfaceEval,
  NurbsInverse,
  TriBounds,
  TriBspline,
  TriCylinder,
  TriSphere,
  TriToroidal,
  TriConical,
  TriRevolution,
  // The advanced-BREP SURFACE tessellator for a surface_of_linear_extrusion
  // face (mesh_utils.h). Not the solid sweep -- that is ExtrudeSolid below.
  // Same word, disjoint call graphs; see geometry-memory-coverage.md.
  TriExtrusion,
  // --- scope kinds: one per instrumented unit of geometry work -------------
  AdvancedFace,
  ExtrudeSolid,
  SweepSolid,
  CsgBoolean,
  // --- sites inside those scopes -------------------------------------------
  ExtrudeCap,
  CsgOperandPrep,
  CsgKernel,
  // The global VertexWelder (Geometry.cpp). Tagged because its containers are
  // clear()/reserve()d and never shrink_to_fit-ed, so capacity it grows inside
  // a scope is still live at scope close and is booked as retained OUTPUT when
  // it is reusable scratch. Which scope kinds that actually reaches is a
  // question this tag answers directly; see geometry-memory-coverage.md.
  VertexWeld,
  Count
};

/** RAII scope marking one unit of instrumented geometry work, bucketed by
 *  `kind` (see the scope-kind block in AllocSite). Nesting is ignored
 *  (outermost scope wins) so helper paths that recurse don't double-count --
 *  which also means a nested scope's kind is discarded, not merged. */
class AllocTelemetryScope {
 public:
  explicit AllocTelemetryScope(AllocSite kind);
  ~AllocTelemetryScope();

  AllocTelemetryScope(const AllocTelemetryScope&) = delete;
  AllocTelemetryScope& operator=(const AllocTelemetryScope&) = delete;

 private:
  bool outermost_ = false;
};

/** RAII: set the active allocation-attribution site for the current thread,
 *  restoring the previous one on exit (nests). Allocations made while active
 *  are counted against this site in the telemetry breakdown. */
class AllocTagScope {
 public:
  explicit AllocTagScope(AllocSite site);
  ~AllocTagScope();

  AllocTagScope(const AllocTagScope&) = delete;
  AllocTagScope& operator=(const AllocTagScope&) = delete;

 private:
  AllocSite previous_;
};

/** Aggregate counters for one scope kind. Read back rather than parsed out of
 *  the stderr report so the ownership and lifetime invariants can be asserted
 *  against known-answer scopes in a unit test (test/alloc_telemetry_test.cpp).
 *
 *  Two exact identities hold over these fields, and both are asserted there:
 *
 *    cumulativeBytes == ownedBytes       + unownedBytes
 *    ownedBytes      == diedBytes        + escapedBytes
 *    allocCalls      == ownedAllocs      + unownedAllocs
 *    unownedAllocs   == unownedFullAllocs + unownedNoTableAllocs
 *    freeCalls       == diedCalls        + foreignFrees
 *
 *  Two counters sit deliberately OUTSIDE all of them, for the same reason:
 *  each records a call that really happened but that neither moved bytes nor
 *  produced anything to classify, so folding either in would inflate a number
 *  the identities rest on.
 *
 *    - `failedAllocs` — an allocation call that returned null.
 *    - `nullFrees` — `free(nullptr)`, which is legal, common, and releases
 *      nothing. `freeCalls` therefore counts non-null frees only, which is
 *      what makes the last identity above exact.
 *
 *  `unownedBytes` is allocation the ownership table could not take. It is
 *  still counted in `allocCalls`, `cumulativeBytes` and the load-wide
 *  denominator, but it cannot be lifetime-classified, so it is a residual
 *  rather than being folded into `died` or `escaped`. Its two causes are
 *  reported apart because their remedies are opposite: `unownedFull*` means
 *  the table hit its capacity and wants a BIGGER
 *  CONWAY_ALLOC_TELEMETRY_TABLE_BITS, while `unownedNoTable*` means the
 *  per-thread table could not be allocated at all and wants a SMALLER one (or
 *  less memory pressure). Counts and the denominator stay valid either way;
 *  only classification is lost. */
struct AllocTelemetryKindTotals {
  uint64_t scopes;
  uint64_t allocCalls;
  uint64_t failedAllocs;
  uint64_t cumulativeBytes;
  uint64_t ownedAllocs;
  uint64_t ownedBytes;
  uint64_t unownedAllocs;
  uint64_t unownedBytes;
  uint64_t unownedFullAllocs;
  uint64_t unownedFullBytes;
  uint64_t unownedNoTableAllocs;
  uint64_t unownedNoTableBytes;
  uint64_t diedCalls;
  uint64_t diedBytes;
  uint64_t maxDiedBytes;
  uint64_t escapedBytes;
  uint64_t maxEscapedBytes;
  uint64_t peakBytesTotal;
  uint64_t maxPeakBytes;
  uint64_t freeCalls;
  uint64_t nullFrees;
  uint64_t freedBytes;
  uint64_t foreignFrees;
  uint64_t foreignBytes;
};

/** Per-(scope kind, site) allocation and death totals; escaped bytes for a
 *  site are `bytes - diedBytes` less that site's share of `unownedBytes`,
 *  which the report prints. */
struct AllocTelemetrySiteTotals {
  uint64_t allocs;
  uint64_t bytes;
  uint64_t diedCalls;
  uint64_t diedBytes;
};

/** Allocator traffic for the whole load, counted outside and inside every
 *  scope alike -- the denominator that turns "this path makes N calls" into a
 *  share. A realloc counts as one free and one allocation. */
struct AllocTelemetryLoadTotals {
  uint64_t allocCalls;
  uint64_t allocFailed;
  uint64_t allocBytes;
  uint64_t freeCalls;
  uint64_t freeNull;
  uint64_t freeBytes;
};

AllocTelemetryKindTotals GetAllocTelemetryKindTotals(AllocSite kind);
AllocTelemetrySiteTotals GetAllocTelemetrySiteTotals(AllocSite kind,
                                                     AllocSite site);
AllocTelemetryLoadTotals GetAllocTelemetryLoadTotals();

/** Print the aggregate histogram/summary to stderr (typically once per model,
 *  from the processor destructor or an explicit binding). */
void DumpAllocTelemetry(const char* label);

/** Reset all aggregate counters (e.g. between models in one process). */
void ResetAllocTelemetry();

#else  // !CONWAY_ALLOC_TELEMETRY — zero-cost stubs

enum class AllocSite {
  Other = 0,
  Earcut,
  Cdt,
  SurfaceEval,
  NurbsInverse,
  TriBounds,
  TriBspline,
  TriCylinder,
  TriSphere,
  TriToroidal,
  TriConical,
  TriRevolution,
  // The advanced-BREP SURFACE tessellator for a surface_of_linear_extrusion
  // face (mesh_utils.h). Not the solid sweep -- that is ExtrudeSolid below.
  // Same word, disjoint call graphs; see geometry-memory-coverage.md.
  TriExtrusion,
  // --- scope kinds: one per instrumented unit of geometry work -------------
  AdvancedFace,
  ExtrudeSolid,
  SweepSolid,
  CsgBoolean,
  // --- sites inside those scopes -------------------------------------------
  ExtrudeCap,
  CsgOperandPrep,
  CsgKernel,
  // The global VertexWelder (Geometry.cpp). Tagged because its containers are
  // clear()/reserve()d and never shrink_to_fit-ed, so capacity it grows inside
  // a scope is still live at scope close and is booked as retained OUTPUT when
  // it is reusable scratch. Which scope kinds that actually reaches is a
  // question this tag answers directly; see geometry-memory-coverage.md.
  VertexWeld,
  Count
};

class AllocTelemetryScope {
 public:
  explicit AllocTelemetryScope(AllocSite) {}
};

class AllocTagScope {
 public:
  explicit AllocTagScope(AllocSite) {}
};

struct AllocTelemetryKindTotals {
  uint64_t scopes;
  uint64_t allocCalls;
  uint64_t failedAllocs;
  uint64_t cumulativeBytes;
  uint64_t ownedAllocs;
  uint64_t ownedBytes;
  uint64_t unownedAllocs;
  uint64_t unownedBytes;
  uint64_t unownedFullAllocs;
  uint64_t unownedFullBytes;
  uint64_t unownedNoTableAllocs;
  uint64_t unownedNoTableBytes;
  uint64_t diedCalls;
  uint64_t diedBytes;
  uint64_t maxDiedBytes;
  uint64_t escapedBytes;
  uint64_t maxEscapedBytes;
  uint64_t peakBytesTotal;
  uint64_t maxPeakBytes;
  uint64_t freeCalls;
  uint64_t nullFrees;
  uint64_t freedBytes;
  uint64_t foreignFrees;
  uint64_t foreignBytes;
};

struct AllocTelemetrySiteTotals {
  uint64_t allocs;
  uint64_t bytes;
  uint64_t diedCalls;
  uint64_t diedBytes;
};

struct AllocTelemetryLoadTotals {
  uint64_t allocCalls;
  uint64_t allocFailed;
  uint64_t allocBytes;
  uint64_t freeCalls;
  uint64_t freeNull;
  uint64_t freeBytes;
};

inline AllocTelemetryKindTotals GetAllocTelemetryKindTotals(AllocSite) {
  return AllocTelemetryKindTotals{};
}

inline AllocTelemetrySiteTotals GetAllocTelemetrySiteTotals(AllocSite,
                                                            AllocSite) {
  return AllocTelemetrySiteTotals{};
}

inline AllocTelemetryLoadTotals GetAllocTelemetryLoadTotals() {
  return AllocTelemetryLoadTotals{};
}

inline void DumpAllocTelemetry(const char*) {}
inline void ResetAllocTelemetry() {}

#endif

}  // namespace conway
