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

  /**
   * Reify this to float geometry and indices for use in a graphics API
   *
   * This is called with a 0 offset automatically by GetVertexData and GetIndexData 
   * and their equivalent sizes if reification has not yet been performed.
   *
   * To force reification with an offset, call this with an offset and it will
   * re-reify if the offset is different than the previously used one.
   *
   * @param offset The offset in space to subtract from each vertex
   */
  reify( offset: Vector3 ): void

  clearReification(): void
}
