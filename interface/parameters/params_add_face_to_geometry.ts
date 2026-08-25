import { Bound3DObject } from '../bound_3D_object'
import { StdVector } from '../std_vector'
import { SurfaceObject } from './surface_object'


/**
 * Parameters for adding a b-rep face to geometry.
 */
export interface ParamsAddFaceToGeometry {
  boundsArray: StdVector<Bound3DObject> // std::vector<IfcBound3D>
  advancedBrep: boolean
  surface: SurfaceObject // IfcSurface
  scaling: number

  /**
   * Diagonal of the extent of the representation that DEFINES this face, in
   * the same units the bound points are in, pinned once per representation.
   * Floors the per-face deflection target at a fraction of itself
   * (`OBJECT_DEFLECTION_FLOOR_FACTOR`), so a face that is one tile of a much
   * larger object is not refined far below a pixel — see
   * bldrs-ai/conway#564. Scoped to the definition rather than the model
   * because tessellation is memoized per representation item; zero leaves
   * the per-face target unfloored.
   */
  representationExtent: number
}
