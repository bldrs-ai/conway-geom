import { Vector3 } from '../vector3'


/** Parameter set for an axis 2 placement 3D */
export interface ParamsAxis2Placement3D {
  position: Vector3
  zAxisRef: Vector3
  xAxisRef: Vector3
  normalizeZ: boolean
  normalizeX: boolean
}
