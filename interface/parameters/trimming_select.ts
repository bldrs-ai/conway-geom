import { Vector2 } from '../vector2'
import { Vector3 } from '../vector3'


export interface TrimmingSelect {
  hasParam: boolean
  hasPos: boolean
  hasLength: boolean
  param: number
  pos: Vector2 | undefined
  pos3D: Vector3 | undefined
}
