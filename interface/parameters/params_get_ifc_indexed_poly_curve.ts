import { StdVector } from '../std_vector'
import { Vector2 } from '../vector2'
import { Vector3 } from '../vector3'
import { Segment } from './segment'


/** Parameters to get an ifc indexed poly curve */
export interface ParamsGetIfcIndexedPolyCurve {
  dimensions: number
  segments: StdVector<Segment>
  points: StdVector<Vector2>
}

export interface ParamsGetIfcIndexedPolyCurve3D {
  dimensions: number
  segments: StdVector<Segment>
  points: StdVector<Vector3>
}
