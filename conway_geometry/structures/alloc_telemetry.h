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
 * What it measures: for every instrumented unit of geometry work (an
 * AllocTelemetryScope), the number of system-allocator calls, the peak live
 * bytes allocated inside the scope, and the bytes still live when it closes.
 * This is exactly the transient footprint a per-thread bump arena must hold,
 * so the aggregate histogram picks the arena size and predicts the spill rate.
 * Counters are thread-local during collection and merged into process-wide
 * atomics on scope exit, so the instrument works identically on the MT build.
 *
 * Two independent axes, and confusing them is the trap this instrument fell
 * into once already (conway#637):
 *
 *   - The SCOPE KIND is which call graph opened the unit -- an advanced-BREP
 *     face, one solid extrusion, one swept solid, one boolean composition.
 *     Every aggregate (scopes, alloc calls, peak, retained, histogram) is
 *     bucketed by it, so paths with different natural units do not blend into
 *     one meaningless average.
 *   - The SITE is which sub-step inside the unit made a given allocation
 *     (AllocTagScope). Sites are attributed within their enclosing scope kind.
 *
 * An allocation is only ever seen when it happens inside an
 * AllocTelemetryScope: the wrappers are inert otherwise. So a call graph with
 * no scope on it reports nothing at all, which reads identically to a call
 * graph that allocates nothing. Until conway#639 audited it, Extrude() and the
 * CSG/boolean path had no scope, and the resulting "zero scoped faces on
 * extrusion models" was read as a measurement of those paths when it was only
 * a statement about where the instrument was. Before believing any null from
 * this instrument, check that a scope actually wraps the path in question.
 */

#include <cstddef>

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

inline void DumpAllocTelemetry(const char*) {}
inline void ResetAllocTelemetry() {}

#endif

}  // namespace conway
