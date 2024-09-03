import { Vector3 } from '../vector3'


/** Parameter set for a 3D cartesian transform operator */
export interface ParamsCartesianTransformationOperator3D {
  position: Vector3
  axis1Ref: Vector3
  axis2Ref: Vector3
  axis3Ref: Vector3
  normalizeAxis1: boolean
  normalizeAxis2: boolean
  normalizeAxis3: boolean
  nonUniform: boolean
  realScale: boolean
  scale1_: number
  scale2_: number
  scale3_: number
}
