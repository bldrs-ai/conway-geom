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
}
