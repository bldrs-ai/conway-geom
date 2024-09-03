import { CurveObject } from './curve_object'


/** Represents an IFC 3D profile object */
export interface IfcProfile3D {
  type: string
  curve: CurveObject
  isConvex: boolean
}
