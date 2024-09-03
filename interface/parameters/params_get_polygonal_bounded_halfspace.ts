import { CurveObject } from '../curve_object'
import { NativeTransform3x3 } from '../native_transform'
import { SurfaceObject } from './surface_object'


/** Parameters for getting a polygonal bounded halfspace */
export interface ParamsGetPolygonalBoundedHalfspace {
  scaleFactor:number
  agreement:boolean
  curve:CurveObject | undefined
  surface:SurfaceObject
  position: NativeTransform3x3 // glm::dmat4
}
