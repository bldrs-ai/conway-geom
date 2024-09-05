import { CurveObject } from '../curve_object'
import { StdVector } from '../std_vector'
import { Vector3 } from '../vector3'


/** Params for getting a loop object */
export interface ParamsGetLoop {
  points: StdVector<Vector3> // std::vector<glm::dvec3>
  edges: StdVector<CurveObject>
}
