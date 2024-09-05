import { NativeTransform3x3 } from '../native_transform'

/** Parameters for getting an ellipse curve */
export interface ParamsGetEllipseCurve {
  radiusX: number
  radiusY: number
  hasPlacement: boolean
  placement: NativeTransform3x3
  circleSegments: number
}
