import { NativeTransform } from '../native_transform'


/** Parameters to get a U shaped curve */
export interface ParamsGetUShapeCurve {
  hasPlacement: boolean
  placement: NativeTransform
  depth: number
  flangeWidth: number
  webThickness: number
  flangeThickness: number
  filletRadius: number
  edgeRadius: number
  flangeScope: number
}
