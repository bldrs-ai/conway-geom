import { Deletable } from './deletable'
import { NativeTransform } from './native_transform'
import { ParseBuffer } from './parse_buffer'
import { Vector3 } from './vector3'


/** A native geometry mesh object */
export interface GeometryObject extends Deletable {
  GetVertexData: () => any
  getPoint(parameter: number): Vector3
  getVertexCount(): number

  normalize(): Vector3
  GetVertexDataSize: () => number
  GetIndexData: () => any
  GetIndexDataSize: () => number
  getAllocationSize(): number
  appendGeometry(parameter: GeometryObject): void
  addComponentTransform(transform: any): void
  appendWithTransform(geometry: GeometryObject, transform: NativeTransform): void
  addComponent(parameter: GeometryObject): void
  clone(): GeometryObject
  applyTransform(parameter: any): void
  dumpToOBJ( preamble: string ): string
  extractVertices( buffer: ParseBuffer ): void
  extractTriangles( buffer: ParseBuffer ): void
  extractVerticesAndTriangles( verticesBuffer: ParseBuffer, trianglesBuffer: ParseBuffer ): void
}
