import { Deletable } from '../deletable'
import { StdVector } from '../std_vector'
import { Vector3 } from '../vector3'


/** The native parameter set for a bspline surface (2D manifold surface in 3D) */
export interface BSplineSurface {
  active: boolean
  uDegree: number
  vDegree: number
  closedU: boolean
  closedV: boolean
  // CurveType has been left out deliberately, it's just metadata -- CS
  // weights have been left out, only weightpoints is set from here -- CS
  controlPoints: StdVector<StdVector<Vector3>>
  uMultiplicity: StdVector<number>
  vMultiplicity: StdVector<number>
  uKnots: StdVector<number>
  vKnots: StdVector<number>
  weightPoints: StdVector<StdVector<number>>
}
