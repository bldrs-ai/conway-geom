import { NativeTransform } from '../native_transform'


/** Parameter set for getting a Z shaped curve */
export interface ParamsGetZShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform
  hasFillet: boolean
  depth: number
  flangeWidth: number
  webThickness: number
  flangeThickness: number
  filletRadius: number
  edgeRadius: number
}
