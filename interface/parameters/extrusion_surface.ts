import { ProfileObject } from '../profile_object'
import { Vector3 } from '../vector3'


/** Parameters for an extrusion surface */
export interface ExtrusionSurface {
  active: boolean
  direction: Vector3
  profile: ProfileObject
  length: number
}
