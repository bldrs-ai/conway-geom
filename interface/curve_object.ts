import { Deletable } from './deletable'
import { OutlineDumpable } from './outline_dumpable'
import { Vector2 } from './vector2'
import { Vector3 } from './vector3'


/**
 * Represents a native curve object that may be 2D or 3D
 *
 * Here a curve is represented explictly as a contiguous set
 * of point samples.
 *
 * A curve may be dumped to wavefront obj or SVG as a string
 * to facilitate debug display.
 */
export interface CurveObject extends Deletable, OutlineDumpable {

  add2d: (coord2D: Vector2) => void
  add3d: (coord3D: Vector3) => void
  getPointsSize: () => number
  get2d: (index: number) => Vector2
  get3d: (index: number) => Vector3
  invert: () => void
  isCCW: () => boolean
  dumpToSVG: ( size: Vector2, offset: Vector2 ) => string
  dumpToOBJ: ( preamble: string ) => string
  clone(): CurveObject
  indices: any
}
