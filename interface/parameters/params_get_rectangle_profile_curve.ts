import { NativeTransform3x3 } from '../native_transform'


/** Parameters for getting a rectangle profile curve */
export interface ParamsGetRectangleProfileCurve  {
  xDim: number
  yDim: number
  hasPlacement: boolean
  matrix: NativeTransform3x3 // glm::dmat3
  thickness: number
}