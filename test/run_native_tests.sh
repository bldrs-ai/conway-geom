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

exit "${STATUS}"
