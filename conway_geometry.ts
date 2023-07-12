
import { default as ConwayGeomWasm } from './Dist/ConwayGeomWasm.js'


export interface GeometryObject {
  getVertexData: () => any
  getVertexDataSize: () => number
  getIndexData: () => any
  getIndexDataSize: () => number
  appendGeometry(parameter: GeometryObject): void
  clone(): GeometryObject
  applyTransform(parameter: any): void
}

export interface CurveObject {
  add2d: () => void
  add3d: () => void
  get2d: (index: number) => any
  get3d: (index: number) => any
  invert: () => void
  isCCW: () => boolean
}

export interface ParamsGetCircleCurve {
  radius: number
  hasPlacement: boolean
  placement: any
}

export interface ParamsCreateNativeIfcProfile {
  curve: CurveObject | undefined
  holes: any | undefined //std::vector<conway::geometry::IfcCurve> 
  isConvex: boolean
  isComposite: boolean
  profiles: any | undefined //std::vector<conway::geometry::IfcProfile>;
}

export interface ParamsGetExtrudedAreaSolid {
  depth: number
  dir: any //glm::dvec3
  profile: any //IfcProfile 
};

export interface ParamsGetAxis2Placement2D {
  isAxis2Placement2D: boolean
  isCartesianTransformationOperator2D: boolean
  isCartesianTransformationOperator2DNonUniform: boolean
  position2D: any
  customAxis1Ref: boolean
  axis1Ref: any
  customAxis2Ref: boolean
  axis2Ref: any
  customScale: boolean
  scale1: number
  customScale2: boolean
  scale2: number
};

export interface Segment {
  isArcType: boolean
  indices: any
}

export interface ParamsGetIfcIndexedPolyCurve {
  dimensions: number
  segments: any
  points: any
};

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
  position: any
  zAxisRef: any
  xAxisRef: any
  normalizeZ: boolean
  normalizeX: boolean
}

/**
 * Internal interface for wasm module, geometry processing
 * OBJ + GLTF + GLB (Draco) Conversions
 */
export class ConwayGeometry {
  public wasmModule?: any
  initialized = false

  /**
   *
   * @param wasmModule_ - Pass loaded wasm module to this function if it's already loaded
   */
  constructor(wasmModule_?: any) {
    if (wasmModule_ !== void 0) {
      this.wasmModule = wasmModule_
    }
  }

  /**
   *
   * @return {Promise<boolean>} - initialization status
   */
  async initialize(): Promise<boolean> {
    if (this.wasmModule === void 0) {
      this.wasmModule = await new ConwayGeomWasm()
    }

    this.initialized = false
    this.initialized = this.wasmModule.initializeGeometryProcessor()

    return this.initialized
  }

  /**
   *
   * @param parameters - ParamsPolygonalFaceSet parsed from data model
   * @return {GeometryObject} - Native geometry object
   */
  getGeometry(parameters: ParamsPolygonalFaceSet): GeometryObject {
    const result = this.wasmModule.getGeometry(parameters)
    return result
  }

  /**
   * 
   * @param parameters - ParamsGetIfcIndexedPolyCurve parsed from data model
   * @returns 
   */
  getIndexedPolyCurve(parameters: ParamsGetIfcIndexedPolyCurve): CurveObject {
    const result = this.wasmModule.getIndexedPolyCurve(parameters)
    return result
  }

  /**
   * 
   * @param parameters ParamsGetCircleCurve parsed from data model
   * @returns 
   */
  getCircleCurve(parameters: ParamsGetCircleCurve): CurveObject {
    const result = this.wasmModule.getCircleCurve(parameters)
    return result
  }

  /**
   * 
   * @param parameters ParamsCreateNativeIfcProfile parsed from data model
   * @returns 
   */
  createNativeIfcProfile(parameters: ParamsCreateNativeIfcProfile): any {
    const result = this.wasmModule.createNativeIfcProfile(parameters)
    return result
  }

  /**
   * 
   * @param parameters ParamsGetExtrudedAreaSolid parsed from data model 
   * @returns 
   */
  getExtrudedAreaSolid(parameters:ParamsGetExtrudedAreaSolid): GeometryObject {
    const result = this.wasmModule.getExtrudedAreaSolid(parameters)
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
    return this.wasmModule.geometryToGltf(geometry, isGlb, outputDraco, fileUri)
  }

  /**
   * 
   * @param parameters - ParamsGetAxis2Placement2D structure
   * @returns {any} - native Axis2Placement2D structure
   */
  getAxis2Placement2D(parameters: ParamsGetAxis2Placement2D): any {
    return this.wasmModule.getAxis2Placement2D(parameters)
  }

  /**
   *
   * @param parameters - ParamsAxis2Placement3D structure
   * @return {any} - native Axis2Placement3D structure
   */
  getAxis2Placement3D(parameters: ParamsAxis2Placement3D): any {
    return this.wasmModule.getAxis2Placement3D(parameters)
  }

  /**
   *
   * @param parameters - ParamsLocalPlacement structure
   * @return {any} = native LocalPlacement structure
   */
  getLocalPlacement(parameters: ParamsLocalPlacement) {
    return this.wasmModule.getLocalPlacement(parameters)
  }

  /**
   *
   * @param geometry - Native Geometry Object
   * @return {string} - containing OBJ file contents
   */
  toObj(geometry: GeometryObject): string {
    return this.wasmModule.geometryToObj(geometry, 0)
  }

  /**
   * Frees the geometry processor
   */
  destroy() {
    this.wasmModule.freeGeometryProcessor()
    this.initialized = false
  }
}

