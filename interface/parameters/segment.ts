import { StdVector } from '../std_vector'


/** Segment in a poly curve */
export interface Segment {
  isArcType: boolean
  indices: StdVector<number>
}
