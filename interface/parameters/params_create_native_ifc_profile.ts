import { CurveObject } from '../curve_object'
import { ProfileObject } from '../profile_object'
import { StdVector } from '../std_vector'


/** Parameters for creating a native IFC profile. */
export interface ParamsCreateNativeIfcProfile {
  curve: CurveObject
  holes?: StdVector< CurveObject > // std::vector<conway::geometry::IfcCurve>
  isConvex: boolean
  isComposite: boolean
  profiles?: StdVector< ProfileObject > // std::vector<conway::geometry::IfcProfile>;
}
