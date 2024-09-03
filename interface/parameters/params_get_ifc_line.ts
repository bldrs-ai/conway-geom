import { Vector2 } from '../vector2'
import { Vector3 } from '../vector3'
import { ParamsGetIfcTrimmedCurve } from './params_get_ifc_trimmed_curve'


/** Parameters get  */
export interface ParamsGetIfcLine {
  dimensions: number
  cartesianPoint2D: Vector2
  cartesianPoint3D: Vector3
  vectorOrientation: Vector3
  vectorMagnitude: number
  isEdge: boolean
  paramsGetIfcTrimmedCurve: ParamsGetIfcTrimmedCurve
}
