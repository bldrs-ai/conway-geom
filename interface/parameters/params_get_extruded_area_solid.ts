import { ProfileObject } from '../profile_object'
import { Vector3 } from '../vector3'


/** Parametrs for getting an extruded area solid */
export interface ParamsGetExtrudedAreaSolid {
  depth: number
  dir: Vector3 // glm::dvec3
  profile: ProfileObject // IfcProfile
}
