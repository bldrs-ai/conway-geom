import { CurveObject } from '../curve_object'


/** Parameters for creating a 3D bound */
export interface ParamsCreateBound3D {
  curve: CurveObject // conway::geometry::IfcCurve
  orientation: boolean
  type: number // uint32_t
}
