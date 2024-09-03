import { NativeTransform3x3 } from '../native_transform'


/** Parameters to get an I shaped curve */
export interface ParamsGetIShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform3x3
  hasFillet: boolean
  width: number
  depth: number
  webThickness: number
  flangeThickness: number
  filletRadius: number
}
