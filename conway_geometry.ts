
import { default as ConwayGeomWasm } from './ConwayGeomWasm.js'


export interface GeometryObject {
    getVertexData: () => any
    getVertexDataSize: () => number
    getIndexData: () => any
    getIndexDataSize: () => number
    addGeometry(parameter: GeometryObject): void
}

export interface ParamsPolygonalFaceSet {
    numPoints: number
    numIndices: number
    indicesPerFace: number
    indexedPolygonalFaceWithVoids: boolean
    points: any
    indices: number[]
}

export interface ResultsGltf {
    success: boolean
    bufferUris: any
    buffers: any
}

/**
 * Internal interface for wasm module, geometry processing
 * OBJ + GLTF + GLB (Draco) Conversions
 */
export class ConwayGeometry {
  modelId: number = -1
  public wasmModule: undefined | any = undefined
  initialized = false

  /**
   *
   * @param wasmModule_ - Pass loaded wasm module to this function if it's already loaded
   */
  constructor(wasmModule_?: any) {
    if (typeof wasmModule_ !== 'undefined') {
      this.wasmModule = wasmModule_
    }
  }

  /**
   *
   * @return {Promise<number>} - modelId after initialization is done
   */
  async initialize(): Promise<number> {
    if (!this.wasmModule) {
      this.wasmModule = await new ConwayGeomWasm()
    }

    this.modelId = this.wasmModule.initializeGeometryProcessor()
    this.initialized = true
    return this.modelId
  }

  /**
   *
   * @param parameters - ParamsPolygonalFaceSet parsed from data model
   * @return {GeometryObject} - Native geometry object
   */
  getGeometry(parameters: ParamsPolygonalFaceSet): GeometryObject {
    const result = this.wasmModule.getGeometry(this.modelId, parameters)
    return result
  }

  /**
   *
   * @param geometry - Native geometry object
   * @param isGlb  - boolean if the output should be a single GLB file
   * @param outputDraco - boolean should the output use Draco compression
   * @param fileUri - string of base name for output files
   * @return {ResultsGltf} - boolean success + buffers + file uris
   */
  toGltf(geometry: GeometryObject, isGlb: boolean, outputDraco: boolean, fileUri: string):
  ResultsGltf {
    return this.wasmModule.geometryToGltf(this.modelId, geometry, isGlb, outputDraco, fileUri)
  }

  /**
   *
   * @param geometry - Native Geometry Object
   * @return {string} - containing OBJ file contents
   */
  toObj(geometry: GeometryObject): string {
    return this.wasmModule.geometryToObj(this.modelId, geometry, 0)
  }

  /**
   * Frees the geometry processor
   */
  destroy() {
    this.wasmModule.freeGeometryProcessor(this.modelId)
    this.modelId = -1
    this.initialized = false
  }
}

