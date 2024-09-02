import { Deletable } from '../deletable'


/**
 * Parameter object for getting a poly curve.
 */
export interface ParamsGetPolyCurve extends Deletable {
  points:number
  pointsLength:number
  dimensions:number
  senseAgreement:boolean
  isEdge:boolean
}
