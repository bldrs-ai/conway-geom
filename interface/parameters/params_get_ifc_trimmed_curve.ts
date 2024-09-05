import { Vector2 } from '../vector2'
import { Vector3 } from '../vector3'


/** Parameter set to get an IFC trimmed curve */
export interface ParamsGetIfcTrimmedCurve {
  masterRepresentation: number
  dimensions: number
  senseAgreement: boolean
  trim1Cartesian2D?: Vector2
  trim1Cartesian3D?: Vector3
  trim1Double: number
  trim2Cartesian2D?: Vector2
  trim2Cartesian3D?: Vector3
  trim2Double: number
  trimExists: boolean
}
