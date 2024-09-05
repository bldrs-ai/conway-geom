import { NativeTransform3x3 } from '../native_transform'


/** Parameters for getting a circle curve */
export interface ParamsGetCircleCurve {
  radius: number
  hasPlacement: boolean
  placement: NativeTransform3x3
  thickness: number
}
