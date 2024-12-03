import { CurveObject } from "../curve_object"

/** Parameter set for getting a swept disk solid */
export interface ParamsGetSweptDiskSolid {
  directrix: CurveObject
  radius: number
  innerRadius: number
  startParam: number
  endParam: number
  closed: boolean
  circleSegments:number
  scalingFactor: number
}
