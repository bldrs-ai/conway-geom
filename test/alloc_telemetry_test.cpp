/*
 * Known-answer tests for the AFTP allocation instrument's ownership, lifetime
 * and load-wide counters (bldrs-ai/conway#653).
 *
 * These are ground truth, not smoke: each scope below allocates a stated
 * number of blocks of stated size and frees a stated subset inside the scope,
 * so every counter has one arithmetically correct value and the test asserts
 * that value rather than "something nonzero appeared". That matters because
 * the failure this instrument keeps producing is a counter that reads zero for
 * a structural reason and gets believed -- two of the five retractions in
 * `design/new/geometry-memory-coverage.md` are that error, and #637's original
 * defect was the same one a level down.
 *
 * The three tests named for retractions are the falsifiable half. Reverting
 * `onFree` to the pre-#653 "subtract every free, clamp at zero" turns each of
 * them red with a specific wrong number, stated in each test's comment, not
 * merely "a failure".
 *
 * WHY THE EXACT COUNTS ARE SOUND HERE. `--wrap` is a link-time rewrite of
 * references in the objects handed to the linker, so only `malloc`/`free`
 * called from THIS translation unit and from alloc_telemetry.cpp are
 * intercepted; libstdc++'s `operator new` and libc's own internal allocations
 * live in shared objects and are not. `testEmptyScopeSeesNoStrayAllocations`
 * asserts that empirically before any of the exact counts below rely on it --
 * if the link ever changes (a static libstdc++, say) that canary fails first
 * and names the reason, rather than every other count drifting by an unknown
 * amount.
 *
 * Standalone by design: includes the instrument's own header, links
 * alloc_telemetry.cpp and nothing else. Built by test/run_native_tests.sh with
 * -DCONWAY_ALLOC_TELEMETRY, -DCONWAY_ALLOC_TELEMETRY_TABLE_BITS=8 and the four
 * --wrap flags; the small table is what makes the overflow behaviour reachable
 * in a test rather than only on a 200 MB model.
 */
#include "conway_geometry/structures/alloc_telemetry.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <malloc.h>

// The overflow test asserts an exact owned/unowned split, which it can only do
// against a known table size. A mismatch here would otherwise surface as an
// arithmetic failure in one test with no indication that the build flag, not
// the instrument, was what moved.
#if !defined( CONWAY_ALLOC_TELEMETRY_TABLE_BITS ) || \
    CONWAY_ALLOC_TELEMETRY_TABLE_BITS != 8
#error "build this test with -DCONWAY_ALLOC_TELEMETRY_TABLE_BITS=8 \
(test/run_native_tests.sh does)"
#endif

using conway::AllocSite;
using conway::AllocTelemetryKindTotals;
using conway::AllocTelemetryLoadTotals;
using conway::AllocTelemetryScope;
using conway::AllocTelemetrySiteTotals;
using conway::AllocTagScope;

namespace {

int failures = 0;

void check( bool condition, const char* what ) {

  if ( condition ) {
    printf( "  ok    %s\n", what );
    return;
  }

  printf( "  FAIL  %s\n", what );
  ++failures;
}

void checkEq( uint64_t actual, uint64_t expected, const char* what ) {

  if ( actual == expected ) {
    printf( "  ok    %s (%llu)\n", what,
            static_cast< unsigned long long >( actual ) );
    return;
  }

  printf( "  FAIL  %s: got %llu, expected %llu\n", what,
          static_cast< unsigned long long >( actual ),
          static_cast< unsigned long long >( expected ) );
  ++failures;
}

// kTableLoadLimit in alloc_telemetry.cpp is (slots / 8) * 7, and the runner
// builds this file with CONWAY_ALLOC_TELEMETRY_TABLE_BITS=8. Both halves of
// that arithmetic are restated here so the overflow test asserts an exact
// split rather than "some allocations were refused".
constexpr uint64_t TABLE_SLOTS = 256;
constexpr uint64_t TABLE_LOAD_LIMIT = ( TABLE_SLOTS / 8 ) * 7;

/**
 * Sum of the allocator's usable sizes for a run of pointers.
 *
 * The instrument books `malloc_usable_size`, not the requested size, because
 * that is the memory the allocator actually took; the expectation has to be
 * built the same way or every byte assertion is off by the rounding.
 *
 * @param pointers Blocks to measure.
 * @param from First index, inclusive.
 * @param to Last index, exclusive.
 * @return {uint64_t} Total usable bytes.
 */
uint64_t usableBytes( void* const* pointers, size_t from, size_t to ) {

  uint64_t total = 0;

  for ( size_t where = from; where < to; ++where ) {
    total += malloc_usable_size( pointers[ where ] );
  }

  return total;
}

/**
 * A scope that allocates nothing must report allocating nothing.
 *
 * This is the canary for every exact count below: it establishes that no
 * allocation from the C++ runtime or from the instrument's own machinery is
 * attributed to an open scope in this binary. The ownership table itself is
 * taken from __real_malloc precisely so it does not land here.
 */
void testEmptyScopeSeesNoStrayAllocations() {

  printf( "\n=== empty scope sees no stray allocations ===\n" );
  conway::ResetAllocTelemetry();

  {
    AllocTelemetryScope scope( AllocSite::AdvancedFace );
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::AdvancedFace );

  checkEq( totals.scopes, 1, "the scope was recorded" );
  checkEq( totals.allocCalls, 0, "no allocation was attributed to it" );
  checkEq( totals.freeCalls, 0, "no free was attributed to it" );
  checkEq( totals.cumulativeBytes, 0, "no bytes were attributed to it" );
}

/**
 * 32 blocks in, 20 freed inside the scope: died is those 20, escaped is the
 * other 12, and the two sum to what the scope allocated.
 *
 * This is the capability #653 asks for in one assertion set. Before it, the
 * escaped figure was whatever `liveBytes` happened to hold after every free in
 * the scope had been subtracted from it regardless of provenance.
 */
void testDiedAndEscapedPartitionTheScope() {

  printf( "\n=== died + escaped partition an owned scope ===\n" );
  conway::ResetAllocTelemetry();

  constexpr size_t BLOCKS = 32;
  constexpr size_t FREED_IN_SCOPE = 20;

  void* blocks[ BLOCKS ] = {};
  uint64_t diedBytes = 0;
  uint64_t escapedBytes = 0;
  uint64_t peakBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::ExtrudeSolid );

    for ( size_t where = 0; where < BLOCKS; ++where ) {
      blocks[ where ] = malloc( 1024 );
    }

    // Measured while every block is still live, so the expectation for the
    // live peak is the whole set -- the frees below all come after.
    peakBytes = usableBytes( blocks, 0, BLOCKS );
    diedBytes = usableBytes( blocks, 0, FREED_IN_SCOPE );
    escapedBytes = usableBytes( blocks, FREED_IN_SCOPE, BLOCKS );

    for ( size_t where = 0; where < FREED_IN_SCOPE; ++where ) {
      free( blocks[ where ] );
      blocks[ where ] = nullptr;
    }
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::ExtrudeSolid );

  checkEq( totals.scopes, 1, "one unit" );
  checkEq( totals.allocCalls, BLOCKS, "every allocation was counted" );
  checkEq( totals.ownedAllocs, BLOCKS, "every allocation was owned" );
  checkEq( totals.unownedAllocs, 0, "none overflowed the table" );
  checkEq( totals.diedCalls, FREED_IN_SCOPE, "in-scope frees were owned" );
  checkEq( totals.foreignFrees, 0, "no free was misread as foreign" );
  checkEq( totals.diedBytes, diedBytes, "died-in-scope bytes" );
  checkEq( totals.escapedBytes, escapedBytes, "escaped bytes" );
  checkEq( totals.maxPeakBytes, peakBytes, "live peak is the full set" );
  checkEq( totals.cumulativeBytes, diedBytes + escapedBytes,
           "cumulative == died + escaped" );
  checkEq( totals.ownedBytes, totals.diedBytes + totals.escapedBytes,
           "owned == died + escaped" );

  // Freeing the survivors outside the scope must not reach back into the
  // closed unit's figures.
  for ( size_t where = FREED_IN_SCOPE; where < BLOCKS; ++where ) {
    free( blocks[ where ] );
  }

  const AllocTelemetryKindTotals after =
    conway::GetAllocTelemetryKindTotals( AllocSite::ExtrudeSolid );

  checkEq( after.escapedBytes, escapedBytes,
           "post-scope frees leave the closed unit alone" );
}

/**
 * A free, inside the scope, of a block allocated before it began, and LARGER
 * than anything the scope holds. Nothing may be subtracted.
 *
 * Red on a revert of ownership: the pre-#653 onFree took this branch's clamp,
 * reporting escapedBytes 0 (instead of the small block's usable size) and one
 * clamped free. This is the mechanism behind the retracted `peak:retained`
 * ratios -- 878,666 of these fired in `csg_boolean` on D3D.
 */
void testLargePreScopeFreeIsNotSubtracted() {

  printf( "\n=== a large pre-scope free is counted, not subtracted ===\n" );
  conway::ResetAllocTelemetry();

  void* preScope = malloc( 64 * 1024 );
  const uint64_t preScopeBytes = malloc_usable_size( preScope );

  void* inScope = nullptr;
  uint64_t inScopeBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::SweepSolid );

    inScope = malloc( 256 );
    inScopeBytes = malloc_usable_size( inScope );
    free( preScope );
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::SweepSolid );

  check( preScopeBytes > inScopeBytes,
         "the pre-scope block is the larger one (the clamping case)" );
  checkEq( totals.allocCalls, 1, "only the in-scope allocation was counted" );
  checkEq( totals.foreignFrees, 1, "the pre-scope free was seen" );
  checkEq( totals.foreignBytes, preScopeBytes, "and measured" );
  checkEq( totals.diedBytes, 0, "it did not count as an in-scope death" );
  checkEq( totals.escapedBytes, inScopeBytes,
           "escaped is the in-scope block, undamaged" );
  checkEq( totals.maxPeakBytes, inScopeBytes, "so is the live peak" );

  free( inScope );
}

/**
 * The same shape with the sizes swapped: a pre-scope block SMALLER than what
 * the scope holds live.
 *
 * This is the case the old clamp could not even see. `size > tls.liveBytes`
 * was false, so the subtraction happened silently and the clamp counter stayed
 * at zero -- which is why "zero clamped frees" was withdrawn as an
 * ownership all-clear. Red on a revert: escapedBytes comes back as
 * inScopeBytes - preScopeBytes with no counter anywhere indicating it.
 */
void testSmallPreScopeFreeIsNotSubtractedSilently() {

  printf( "\n=== a small pre-scope free is not subtracted silently ===\n" );
  conway::ResetAllocTelemetry();

  void* preScope = malloc( 256 );
  const uint64_t preScopeBytes = malloc_usable_size( preScope );

  void* inScope = nullptr;
  uint64_t inScopeBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::CsgBoolean );

    inScope = malloc( 64 * 1024 );
    inScopeBytes = malloc_usable_size( inScope );
    free( preScope );
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::CsgBoolean );

  check( preScopeBytes < inScopeBytes,
         "the pre-scope block is the smaller one (the silent case)" );
  checkEq( totals.foreignFrees, 1, "the pre-scope free was seen" );
  checkEq( totals.foreignBytes, preScopeBytes, "and measured" );
  checkEq( totals.escapedBytes, inScopeBytes,
           "escaped is the in-scope block, not it minus the foreign free" );

  free( inScope );
}

/**
 * Allocations made with no scope open reach the load-wide counters and no
 * others.
 *
 * The counter this pins is the one that turns "csg_boolean makes 43.2 M
 * allocator calls" into a share. Its validation has to be that it counts what
 * the per-scope instrument does NOT: asserting it merely moved would pass with
 * a counter that only fires inside scopes, which is exactly the instrument
 * conway#637 already had.
 */
void testLoadWideDenominatorCountsUnscopedTraffic() {

  printf( "\n=== the load-wide denominator counts unscoped traffic ===\n" );
  conway::ResetAllocTelemetry();

  constexpr size_t BLOCKS = 64;

  void* blocks[ BLOCKS ] = {};

  for ( size_t where = 0; where < BLOCKS; ++where ) {
    blocks[ where ] = malloc( 512 );
  }

  const uint64_t allocated = usableBytes( blocks, 0, BLOCKS );

  for ( size_t where = 0; where < BLOCKS; ++where ) {
    free( blocks[ where ] );
  }

  const AllocTelemetryLoadTotals load = conway::GetAllocTelemetryLoadTotals();

  checkEq( load.allocCalls, BLOCKS, "unscoped allocations were counted" );
  checkEq( load.allocBytes, allocated, "unscoped bytes were counted" );
  checkEq( load.freeCalls, BLOCKS, "unscoped frees were counted" );
  checkEq( load.freeBytes, allocated, "unscoped freed bytes were counted" );

  uint64_t scopedCalls = 0;

  for ( int kind = 0; kind < static_cast< int >( AllocSite::Count ); ++kind ) {
    scopedCalls += conway::GetAllocTelemetryKindTotals(
        static_cast< AllocSite >( kind ) ).allocCalls;
  }

  checkEq( scopedCalls, 0,
           "and none of it was attributed to any scope kind" );

  // The denominator has to remain a superset once a scope IS open, or a share
  // taken against it could exceed 100 %.
  //
  // The usable size is read inside the scope rather than after it for a
  // reason: gcc elides an allocate/free pair whose pointer is otherwise
  // unused, and it did elide this one, which reported the denominator failing
  // to count an in-scope allocation that had never happened. Reading the size
  // makes the block observable, and the assertion stronger.
  void* scoped = nullptr;
  uint64_t scopedBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::ExtrudeSolid );

    scoped = malloc( 4096 );
    scopedBytes = malloc_usable_size( scoped );
  }

  const AllocTelemetryLoadTotals afterScope =
    conway::GetAllocTelemetryLoadTotals();
  const AllocTelemetryKindTotals kindTotals =
    conway::GetAllocTelemetryKindTotals( AllocSite::ExtrudeSolid );

  checkEq( afterScope.allocCalls, load.allocCalls + 1,
           "an in-scope allocation also reaches the denominator" );
  checkEq( afterScope.allocBytes, load.allocBytes + scopedBytes,
           "with its bytes" );
  checkEq( kindTotals.allocCalls, 1, "and it reached the scope kind too" );
  check( afterScope.allocCalls > kindTotals.allocCalls,
         "the denominator stays a strict superset of in-scope calls" );

  free( scoped );
}

/**
 * More live in-scope pointers than the ownership table can hold.
 *
 * The stated policy is counted-but-unowned, never silent: the refused
 * allocations still appear in allocCalls, in cumulativeBytes and in the
 * load-wide denominator, and the residual is carried as its own field so the
 * report can say the measurement degraded instead of folding it into died or
 * escaped.
 */
void testOverflowIsCountedNotSilent() {

  printf( "\n=== a full ownership table degrades loudly ===\n" );
  conway::ResetAllocTelemetry();

  constexpr size_t BLOCKS = 300;

  void* blocks[ BLOCKS ] = {};

  {
    AllocTelemetryScope scope( AllocSite::CsgBoolean );

    for ( size_t where = 0; where < BLOCKS; ++where ) {
      blocks[ where ] = malloc( 64 );
    }
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::CsgBoolean );
  const AllocTelemetryLoadTotals load = conway::GetAllocTelemetryLoadTotals();
  const uint64_t allocated = usableBytes( blocks, 0, BLOCKS );

  check( BLOCKS > TABLE_LOAD_LIMIT, "the run really does overflow the table" );
  checkEq( totals.allocCalls, BLOCKS, "every allocation was still counted" );
  checkEq( totals.ownedAllocs, TABLE_LOAD_LIMIT, "ownership stopped at the limit" );
  checkEq( totals.unownedAllocs, BLOCKS - TABLE_LOAD_LIMIT,
           "the remainder is reported as unowned" );
  checkEq( totals.cumulativeBytes, allocated,
           "cumulative volume is unaffected by the table" );
  checkEq( totals.ownedBytes + totals.unownedBytes, totals.cumulativeBytes,
           "cumulative == owned + unowned" );
  checkEq( load.allocCalls, BLOCKS,
           "and the load-wide denominator is unaffected too" );

  for ( size_t where = 0; where < BLOCKS; ++where ) {
    free( blocks[ where ] );
  }
}

/**
 * Sustained allocate/free churn with a stable live set.
 *
 * This is the test for the table's deletion path. `tableErase` compacts the
 * probe chain behind a removed entry rather than leaving a tombstone, which is
 * what keeps occupancy proportional to LIVE allocations rather than to total
 * allocations -- a `csg_boolean` composition makes about 6,700 allocations and
 * holds far fewer at once. If the compaction moved an entry past its own home
 * slot the later lookup would miss, and the miss would surface as a foreign
 * free on a pointer the scope certainly owns. So foreignFrees == 0 over 5,000
 * churn cycles is the assertion that the deletion is correct, and the exact
 * died/escaped split is the assertion it stayed correct.
 */
void testChurnKeepsOwnershipExact() {

  printf( "\n=== churn keeps ownership exact ===\n" );
  conway::ResetAllocTelemetry();

  constexpr size_t LIVE = 150;
  constexpr size_t CYCLES = 5000;
  constexpr size_t FREED_AT_END = 100;

  void* live[ LIVE ] = {};
  uint64_t diedBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::CsgBoolean );

    for ( size_t where = 0; where < LIVE; ++where ) {
      live[ where ] = malloc( 32 + ( where % 7 ) * 16 );
    }

    // A plain LCG: the point is a non-monotonic eviction order, which is what
    // produces the long probe chains compaction has to repair, and it must not
    // itself allocate.
    uint32_t state = 12345;

    for ( size_t cycle = 0; cycle < CYCLES; ++cycle ) {
      state = state * 1664525u + 1013904223u;

      const size_t slot = state % LIVE;

      diedBytes += malloc_usable_size( live[ slot ] );
      free( live[ slot ] );
      live[ slot ] = malloc( 32 + ( cycle % 11 ) * 16 );
    }

    for ( size_t where = 0; where < FREED_AT_END; ++where ) {
      diedBytes += malloc_usable_size( live[ where ] );
      free( live[ where ] );
      live[ where ] = nullptr;
    }
  }

  const AllocTelemetryKindTotals totals =
    conway::GetAllocTelemetryKindTotals( AllocSite::CsgBoolean );

  checkEq( totals.allocCalls, LIVE + CYCLES, "every allocation was counted" );
  checkEq( totals.ownedAllocs, LIVE + CYCLES,
           "compaction kept occupancy under the limit throughout" );
  checkEq( totals.unownedAllocs, 0, "so nothing overflowed" );
  checkEq( totals.foreignFrees, 0,
           "no owned pointer was lost by the compaction" );
  checkEq( totals.diedCalls, CYCLES + FREED_AT_END, "in-scope deaths" );
  checkEq( totals.diedBytes, diedBytes, "in-scope died bytes" );
  checkEq( totals.escapedBytes, totals.cumulativeBytes - diedBytes,
           "escaped is the remainder" );

  for ( size_t where = FREED_AT_END; where < LIVE; ++where ) {
    free( live[ where ] );
  }
}

/**
 * Lifetime is classified per site as well as per scope kind.
 *
 * This is the half that answers the third retraction: the global VertexWelder
 * grows capacity inside a scope and never shrinks it, so its bytes are escaped
 * but they are reusable scratch rather than that unit's output. Only a
 * per-site split can separate the two, so the split is pinned here.
 */
void testSitesAreLifetimeClassifiedSeparately() {

  printf( "\n=== sites are lifetime-classified separately ===\n" );
  conway::ResetAllocTelemetry();

  void* welded = nullptr;
  uint64_t weldBytes = 0;
  uint64_t kernelBytes = 0;

  {
    AllocTelemetryScope scope( AllocSite::CsgBoolean );

    {
      AllocTagScope tag( AllocSite::CsgKernel );

      void* transient = malloc( 2048 );

      kernelBytes = malloc_usable_size( transient );
      free( transient );
    }

    {
      AllocTagScope tag( AllocSite::VertexWeld );

      welded = malloc( 4096 );
      weldBytes = malloc_usable_size( welded );
    }
  }

  const AllocTelemetrySiteTotals kernel =
    conway::GetAllocTelemetrySiteTotals( AllocSite::CsgBoolean,
                                         AllocSite::CsgKernel );
  const AllocTelemetrySiteTotals weld =
    conway::GetAllocTelemetrySiteTotals( AllocSite::CsgBoolean,
                                         AllocSite::VertexWeld );

  checkEq( kernel.bytes, kernelBytes, "kernel bytes attributed to the kernel" );
  checkEq( kernel.diedBytes, kernelBytes, "and they died in scope" );
  checkEq( weld.bytes, weldBytes, "weld bytes attributed to the weld" );
  checkEq( weld.diedBytes, 0, "and they escaped it" );
  checkEq( weld.bytes - weld.diedBytes, weldBytes,
           "so the persistent cache is separable from the unit's output" );

  free( welded );
}

}  // namespace

int main() {

  testEmptyScopeSeesNoStrayAllocations();
  testDiedAndEscapedPartitionTheScope();
  testLargePreScopeFreeIsNotSubtracted();
  testSmallPreScopeFreeIsNotSubtractedSilently();
  testLoadWideDenominatorCountsUnscopedTraffic();
  testOverflowIsCountedNotSilent();
  testChurnKeepsOwnershipExact();
  testSitesAreLifetimeClassifiedSeparately();

  // Not an assertion -- the report is what a profiling run reads, and running
  // it here is what stops a format-string or divide-by-zero defect in it from
  // first appearing on a 35-second wasm load. The numbers below are the last
  // test's, so they are small and checkable by eye.
  printf( "\n=== report (stderr), from the last test's counters ===\n" );
  fflush( stdout );
  conway::DumpAllocTelemetry( "native ground-truth run" );

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
