import { Deletable } from '../deletable'
import { StdVector } from '../std_vector'


/** Indexed polygonal face */
export interface IndexedPolygonalFace extends Deletable {
  indices: StdVector< number >
  face_starts: StdVector< number >
}
