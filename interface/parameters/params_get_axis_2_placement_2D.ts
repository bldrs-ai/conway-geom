import { Vector2 } from '../vector2'


/** Parameters for getting a native transform from an Axis 2 placement 2D */
export interface ParamsGetAxis2Placement2D {
  isAxis2Placement2D: boolean
  isCartesianTransformationOperator2D: boolean
  isCartesianTransformationOperator2DNonUniform: boolean
  position2D: Vector2
  customAxis1Ref: boolean
  axis1Ref: Vector2
  customAxis2Ref: boolean
  axis2Ref: Vector2
  customScale: boolean
  scale1: number
  customScale2: boolean
  scale2: number
}
