import { NativeTransform3x3 } from '../native_transform'


/** Parameters to get a T shaped curve */
export interface ParamsGetTShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform3x3
  hasFillet: boolean
  depth: number
  width: number
  webThickness: number
  filletRadius: number
  flangeEdgeRadius: number
  // webEdgeRadius:number
  // webSlope:number
  flangeScope: number
}
