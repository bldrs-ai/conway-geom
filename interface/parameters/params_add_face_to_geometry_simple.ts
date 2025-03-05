import { Bound3DObject } from '../bound_3D_object'
import { StdVector } from '../std_vector'


/**
 * Parameters for adding a b-rep face to geometry.
 */
export interface ParamsAddFaceToGeometrySimple {
  boundsArray: StdVector<Bound3DObject> // std::vector<IfcBound3D>
  scaling: number
}
