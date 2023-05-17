
import { default as ConwayGeomWasm } from './Dist/ConwayGeomWasm.js'


export interface GeometryObject {
  getVertexData: () => any
  getVertexDataSize: () => number
  getIndexData: () => any
  getIndexDataSize: () => number
  addGeometry(parameter: GeometryObject): void
  applyTransform(parameter: any): void
}

export interface IndexedPolygonalFace {
  indices: any
  face_starts: any
}

export interface ParamsPolygonalFaceSet {
  indicesPerFace: number
  points: any
  faces: any
}

export interface ResultsGltf {
  success: boolean
  bufferUris: any
  buffers: any
}

export interface ParamsLocalPlacement {
  useRelPlacement: boolean
  axis2Placement: any
  relPlacement: any
}

export interface ParamsAxis2Placement3D {
  position:any;
  zAxisRef:any;
  xAxisRef:any;
  normalizeZ:boolean;
  normalizeX:boolean;
}

// Graph class
class IfcDag {
  adjacencyList: Map<number, number[]>

  constructor() {
    this.adjacencyList = new Map()
  }

  addEdge(product: number, representation: number) {
    if (!this.adjacencyList.has(product)) {
      this.adjacencyList.set(product, [])
    }
    this.adjacencyList.get(product)!.push(representation)
  }
}

/**
 * Internal interface for wasm module, geometry processing
 * OBJ + GLTF + GLB (Draco) Conversions
 */
export class ConwayGeometry {
  modelId: number = -1
  public wasmModule: undefined | any = undefined
  initialized = false
  //map localID of transformation to localID of 1 or more geometries 
  public transformMapping = new Map<number, number[]>()

  public graph: IfcDag = new IfcDag()

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

  getAxis2Placement3D(parameters:ParamsAxis2Placement3D) {
    return this.wasmModule.getAxis2Placement3D(this.modelId, parameters)
  }

  getLocalPlacement(parameters:ParamsLocalPlacement) {
    return this.wasmModule.getLocalPlacement(this.modelId, parameters)
  }

  topologicalSort(graph: IfcDag): number[] | null {
    const result: number[] = []
    const visited = new Map<number, boolean>()
    const visiting = new Map<number, boolean>()

    for (const [node] of graph.adjacencyList) {
      if (!visited.has(node) && !dfs(node)) {
        return null // not a DAG
      }
    }

    function dfs(node: number): boolean {
      visiting.set(node, true)

      for (const neighbor of graph.adjacencyList.get(node) || []) {
        if (visited.has(neighbor)) {
          continue
        }
        if (visiting.has(neighbor) || !dfs(neighbor)) {
          return false
        }
      }

      visiting.delete(node)
      visited.set(node, true)
      result.push(node)

      return true
    }

    return result.reverse()
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

