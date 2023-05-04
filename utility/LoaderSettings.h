/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

namespace webifc::utility {

struct LoaderSettings {
  bool COORDINATE_TO_ORIGIN = false;
  bool USE_FAST_BOOLS = true;  // TODO: This needs to be fixed in the future to
                               // rely on elalish/manifold
  int CIRCLE_SEGMENTS_LOW = 5;
  int CIRCLE_SEGMENTS_MEDIUM = 8;
  int CIRCLE_SEGMENTS_HIGH = 12;
  int BOOL_ABORT_THRESHOLD = 10000;  // 10k verts
  uint32_t TAPE_SIZE = 67108864;     // probably no need for anyone other than
                                     // web-ifc devs to change this
  uint32_t MEMORY_LIMIT = 3221225472;
};
}  // namespace webifc::utility
