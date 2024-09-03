import { NativeTransform3x3 } from '../native_transform'


/** Parameters for getting an L shaped curve */
export interface ParamsGetLShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform3x3
  hasFillet: boolean
  filletRadius: number
  depth: number
  width: number
  thickness: number
  edgeRadius: number
  legSlope: number
  // centerOfGravityInX:number
  // centerOfGravityInY:number
}
