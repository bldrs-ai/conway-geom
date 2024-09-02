import { Deletable } from './deletable'
import { GeometryObject } from './geometry_object'
import { NativeTransform4x4 } from './native_transform'


/**
 * A native collection of geometry objects with associated transforms.
 */
export interface GeometryCollection extends Deletable {

  addComponentWithTransform(geometry: GeometryObject, transform: NativeTransform4x4): void

  materialIndex: number
  hasDefaultMaterial: boolean

  readonly currentSize: number
}
