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
   * Diagonal of the whole model's extent, in the same units the bound points
   * are in, pinned once per load. Floors the per-face deflection target at a
   * fraction of itself (`MODEL_DEFLECTION_FLOOR_FACTOR`), so a face that is
   * one tile of a much larger object is not refined far below a pixel — see
   * bldrs-ai/conway#564. Zero leaves the per-face target unfloored.
   */
  modelExtent: number
}
