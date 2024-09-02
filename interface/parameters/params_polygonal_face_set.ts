import { StdVector } from '../std_vector'
import { Vector3 } from '../vector3'
import { IndexedPolygonalFace } from './indexed_polygonal_face'


/** Parameter set to get a polygonal face set */
export interface ParamsPolygonalFaceSet {
  indicesPerFace: number
  points: StdVector< Vector3 >
  faces: StdVector< IndexedPolygonalFace >
}
