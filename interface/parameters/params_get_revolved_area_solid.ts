import { NativeTransform4x4 } from '../native_transform'
import { ProfileObject } from '../profile_object'
import { Vector3 } from '../vector3'


/** Parametrs for getting a revolved area solid */
export interface ParamsGetRevolvedAreaSolid {
  placement: NativeTransform4x4 // glm::dmat4
  axis: Vector3 // glm::dvec3
  axisPosition: Vector3 // glm::dvec3
  angle:number
  profile: ProfileObject // IfcProfile
  scalingFactor:number
  circleSegments:number
}