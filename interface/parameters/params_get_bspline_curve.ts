import { StdVector } from '../std_vector'
import { Vector2 } from '../vector2'
import { Vector3 } from '../vector3'


/** Parameters for getting a bspline curve */
export interface ParamsGetBSplineCurve {
  dimensions: number
  degree: number
  points2: StdVector<Vector2>
  points3: StdVector<Vector3>
  knots: StdVector<number>
  weights: StdVector<number>
  senseAgreement:boolean
  isEdge:boolean
}
