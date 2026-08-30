#!/usr/bin/env bash
#
# Build and run the header-only native tests.
#
# Deliberately NOT the `conway_geom_native_tests` GENie target: that target does
# not currently build (its include paths are missing utility/), and repairing it
# is not in the scope of the change these tests came with. The tests here are
# self-contained translation units that include the header under test and link
# nothing, so a plain compiler and the repo's own include tree is all they need.
#
# Usage: test/run_native_tests.sh [compiler]

set -euo pipefail

cd "$( dirname "$0" )/.."

COMPILER="${1:-${CXX:-g++}}"
OUT="$( mktemp -d )"
trap 'rm -rf "${OUT}"' EXIT

# The same include tree the wasm build uses, minus the parts only the
# emscripten targets need.
INCLUDES=(
  -I. -Iinclude -Iutility -Iconway_geometry -Ilogging
  -Iexternal/tinynurbs/include
  -Iexternal/manifold/src
  -Iexternal/manifold/src/utilities/include
  -Iexternal/glm
  -Iexternal/earcut.hpp/include
  -Iexternal/manifold/src/collider/include
  -Iexternal/manifold/src/third_party/thrust
  -Iexternal/manifold/src/manifold/include
  -Iexternal/manifold/src/polygon/include
  -Iexternal/manifold/src/sdf/include
  -Iexternal/manifold/src/third_party/graphlite/include
  -Iexternal/manifold/src/third_party/glm
  -Iexternal/CDT/CDT/include
  -Iexternal/fast_float/include
)

STATUS=0

for source in test/*_test.cpp; do

  # Only the self-contained ones; the GENie-target tests need linking.
  case "${source}" in
    test/nurbs_seam_test.cpp | test/ribbon_loft_test.cpp | \
    test/spherical_trim_test.cpp | test/inverse_wrong_sheet_test.cpp ) ;;
    * ) continue ;;
  esac

  name="$( basename "${source}" .cpp )"

  echo "=== ${name} ==="

  "${COMPILER}" -std=c++20 -O1 -DREAL_T_IS_DOUBLE "${INCLUDES[@]}" \
    "${source}" -o "${OUT}/${name}"

  if ! "${OUT}/${name}"; then
    STATUS=1
  fi
done

# The allocation instrument is the one test here that cannot be a
# link-nothing translation unit: it is the --wrap hooks that are under test, so
# it needs alloc_telemetry.cpp compiled alongside it, the compile gate defined,
# and the same four --wrap flags genie.lua adds for the opt-in wasm build.
# GNU ld honours --wrap identically to wasm-ld, so this runs natively.
#
# TABLE_BITS=8 shrinks the per-thread ownership table from 2^18 slots to 2^8,
# which is what puts its overflow behaviour inside reach of a unit test instead
# of only a 200 MB model. The test restates that number and the 7/8 load limit
# to assert an exact owned/unowned split, so the two must be changed together.
echo "=== alloc_telemetry_test ==="

"${COMPILER}" -std=c++20 -O1 -DREAL_T_IS_DOUBLE \
  -DCONWAY_ALLOC_TELEMETRY -DCONWAY_ALLOC_TELEMETRY_TABLE_BITS=8 \
  "${INCLUDES[@]}" \
  test/alloc_telemetry_test.cpp \
  conway_geometry/structures/alloc_telemetry.cpp \
  -o "${OUT}/alloc_telemetry_test" \
  -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=free

if ! "${OUT}/alloc_telemetry_test"; then
  STATUS=1
fi

exit "${STATUS}"
