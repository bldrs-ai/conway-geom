import { Vector2 } from './vector2'
import { WavefrontDumpable } from './wavefront_dumpable'


export interface OutlineDumpable extends WavefrontDumpable {
  dumpToSVG: ( size: Vector2, offset: Vector2 ) => string
}
