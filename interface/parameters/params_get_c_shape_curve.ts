import { NativeTransform3x3 } from '../native_transform'


/** Parameters for getting a C shape curve */
export interface ParamsGetCShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform3x3
  hasFillet: boolean
  depth: number
  width: number
  thickness: number
  girth: number
  filletRadius: number
}
