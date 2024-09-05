import { Deletable } from '../deletable'


/** Parameter set for getting triangulated face set geometry */
export interface ParamsGetTriangulatedFaceSetGeometry extends Deletable {
  points:number /* WASM pointer */
  pointsArrayLength:number
  indices:number /* WASM pointer */
  indicesArrayLength: number
}
