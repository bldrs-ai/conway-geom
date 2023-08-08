
import { default as ConwayGeomWasm } from './Dist/ConwayGeomWasm.js'


export interface StdVector<T> {
  delete(): unknown

  resize(size: number, value?: T): void

  push_back(value: T): void

  size(): number

  empty(): boolean

  set(index: number, value?: T): void
  get(index: number): T
}

export interface GeometryObject {
  getVertexData: () => any
  getVertexDataSize: () => number
  getIndexData: () => any
  getIndexDataSize: () => number
  appendGeometry(parameter: GeometryObject): void
  clone(): GeometryObject
  applyTransform(parameter: any): void
  materialIndex: number
  hasDefaultMaterial: boolean
}

export interface SurfaceObject {

}

export interface ParamsAddFaceToGeometry {
  boundsArray: StdVector<Bound3DObject> // std::vector<IfcBound3D>
  advancedBrep: boolean
  surface: SurfaceObject // IfcSurface
}

export interface ParamsGetLoop {
  isEdgeLoop: boolean
  points: StdVector<Vector3> // std::vector<glm::dvec3>
}

export interface ParamsCreateBound3D {
  curve: CurveObject // conway::geometry::IfcCurve
  orientation: boolean
  type: number // uint32_t
}

export interface Vector4 {
  x: number
  y: number
  z: number
  w: number
}

/* eslint-disable no-shadow,no-unused-vars,no-magic-numbers */
export enum BlendMode {
  OPAQUE = 0,
  BLEND = 1,
  MASK = 2,
}

/**
 * @return {any}
 */
export function toAlphaMode(wasmModule: any, blendMode: BlendMode): any {

  switch (blendMode) {
    case BlendMode.OPAQUE:

      return wasmModule.BlendMode.OPAQUE

    case BlendMode.BLEND:

      return wasmModule.BlendMode.BLEND

    case BlendMode.MASK:

      return wasmModule.BlendMode.MASK

    default: console.log('Wrong argument passed to toAlphaMode()')
  }
}

export interface MaterialObject {

  base: Vector4
  metallic: number
  roughness: number
  alphaCutoff: number

  ior: number

  specular?: Vector4

  alphaMode: any

  doubleSided: boolean
}

export interface Bound3DObject {

}

export interface CurveObject {
  add2d: (coord2D: any) => void
  add3d: (coord3D: any) => void
  getPointsSize: () => number
  get2d: (index: number) => any
  get3d: (index: number) => any
  invert: () => void
  isCCW: () => boolean
}

export interface ProfileObject {
  getType: () => string
  getCurve: () => CurveObject
  getHoles: () => any// CurveObject[];
  isConvex: () => boolean
  isComposite: () => boolean
  getProfiles: () => any// ProfileObject[];
}

export interface ParamsGetCircleCurve {
  radius: number
  hasPlacement: boolean
  placement: any
}

export interface ParamsCreateNativeIfcProfile {
  curve: CurveObject | undefined
  holes: any | undefined // std::vector<conway::geometry::IfcCurve>
  isConvex: boolean
  isComposite: boolean
  profiles: any | undefined // std::vector<conway::geometry::IfcProfile>;
}

export interface ParamsGetExtrudedAreaSolid {
  depth: number
  dir: any // glm::dvec3
  profile: ProfileObject // IfcProfile
}

export interface ParamsGetRectangleProfileCurve {
  xDim: number
  yDim: number
  hasPlacement: boolean
  matrix: any // glm::dmat3
}

export interface ParamsGetHalfspaceSolid {
  flipWinding: boolean,
  optionalLinearScalingFactor: number,
}

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
}

export interface Segment {
  isArcType: boolean
  indices: any
}

export interface ParamsGetIfcIndexedPolyCurve {
  dimensions: number
  segments: any
  points: any
}

export interface ParamsGetIfcCircle {
  dimensions: number
  axis2Placement2D: any
  axis2Placement3D: any
  radius: number
  paramsGetIfcTrimmedCurve: ParamsGetIfcTrimmedCurve
}

export interface ParamsGetIfcTrimmedCurve {
  masterRepresentation: number
  dimensions: number
  senseAgreement: boolean
  trim1Cartesian2D: any
  trim1Cartesian3D: any
  trim1Double: number
  trim2Cartesian2D: any
  trim2Cartesian3D: any
  trim2Double: number
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

export interface Vector3 {
  x: number
  y: number
  z: number
}

export interface Vector2 {
  x: number
  y: number
}

export interface ParamsCartesianTransformationOperator3D {
  position: Vector3
  axis1Ref: Vector3
  axis2Ref: Vector3
  axis3Ref: Vector3
  normalizeAxis1: boolean
  normalizeAxis2: boolean
  normalizeAxis3: boolean
  nonUniform: boolean
  realScale: boolean
  scale1_: number
  scale2_: number
  scale3_: number
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

export interface ParamsGetBooleanResult {
  flatFirstMesh: any
  flatSecondMesh: any
  operatorType: number
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
   * @param parameters ParamsGetLoop parsed from data model
   * @return {CurveObject}
   */
  getLoop(parameters: ParamsGetLoop): CurveObject {
    const result = this.wasmModule.getLoop(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsAddFaceToGeometry parsed from data model
   */
  addFaceToGeometry(parameters: ParamsAddFaceToGeometry, geometry: GeometryObject): void {
    this.wasmModule.addFaceToGeometry(parameters, geometry)
  }

  /**
   *
   * @param parameters - ParamsCartesianTransformationOperator3D parsed from data model
   * @return {GeometryObject} - Native geometry object
   */
  getCartesianTransformationOperator3D(parameters: ParamsCartesianTransformationOperator3D): any {
    const result = this.wasmModule.getCartesianTransformationOperator3D(parameters)
    return result
  }

  /**
   *
   * @param parameters - ParamsPolygonalFaceSet parsed from data model
   * @return {GeometryObject} - Native geometry object
   */
  getPolygonalFaceSetGeometry(parameters: ParamsPolygonalFaceSet): GeometryObject {
    const result = this.wasmModule.getPolygonalFaceSetGeometry(parameters)
    return result
  }

  /**
   *
   * @param parameters
   * @return {CurveObject} - Native Curve Object
   */
  getIfcCircle(parameters: ParamsGetIfcCircle): CurveObject {
    const result = this.wasmModule.getIfcCircle(parameters)
    return result
  }

  /**
   *
   * @param parameters - ParamsGetIfcIndexedPolyCurve parsed from data model
   * @return {CurveObject}
   */
  getIndexedPolyCurve(parameters: ParamsGetIfcIndexedPolyCurve): CurveObject {
    const result = this.wasmModule.getIndexedPolyCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetCircleCurve parsed from data model
   * @return {CurveObject}
   */
  getCircleCurve(parameters: ParamsGetCircleCurve): CurveObject {
    const result = this.wasmModule.getCircleCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsCreateNativeIfcProfile parsed from data model
   * @return {ProfileObject}
   */
  createNativeIfcProfile(parameters: ParamsCreateNativeIfcProfile): ProfileObject {
    const result: ProfileObject = this.wasmModule.createNativeIfcProfile(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsCreateBound3D parsed from data model
   * @return {Bound3DObject}
   */
  createBound3D(parameters: ParamsCreateBound3D): Bound3DObject {
    const result = this.wasmModule.createBound3D(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetHalfspaceSolid parsed from data model
   * @return {GeometryObject}
   */
  getHalfSpaceSolid(parameters: ParamsGetHalfspaceSolid): GeometryObject {
    const result = this.wasmModule.getHalfSpaceSolid(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetRectangleProfileCurve parsed from data model
   * @return {CurveObject}
   */
  getRectangleProfileCurve(parameters: ParamsGetRectangleProfileCurve): CurveObject {
    const result = this.wasmModule.getRectangleProfileCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetExtrudedAreaSolid parsed from data model
   * @return {GeometryObject}
   */
  getExtrudedAreaSolid(parameters: ParamsGetExtrudedAreaSolid): GeometryObject {
    const result = this.wasmModule.getExtrudedAreaSolid(parameters)
    return result
  }

  /**
   *
   * @param parameters
   * @return {GeometryObject}
   */
  getBooleanResult(parameters: ParamsGetBooleanResult): GeometryObject {
    const result = this.wasmModule.getBooleanResult(parameters)
    return result
  }

  /**
   *
   * @param geometry - Vector of native geometry object
   * @param isGlb  - boolean if the output should be a single GLB file
   * @param outputDraco - boolean should the output use Draco compression
   * @param fileUri - string of base name for output files
   * @return {ResultsGltf} - boolean success + buffers + file uris
   */
  toGltf(
      geometry: StdVector<GeometryObject>,
      materials: StdVector<MaterialObject>,
      isGlb: boolean,
      outputDraco: boolean,
      fileUri: string):
    ResultsGltf {
    return this.wasmModule.geometryToGltf(geometry, materials, isGlb, outputDraco, fileUri)
  }

  /**
   *
   * @param parameters - ParamsGetAxis2Placement2D structure
   * @return {any} - native Axis2Placement2D structure
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

