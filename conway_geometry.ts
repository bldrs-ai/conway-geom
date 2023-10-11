
import { LocateFileHandlerFn } from '../../src/shim/ifc_api.js'
import { CanonicalMaterial } from '../../src/core/canonical_material.js'
import { default as ConwayGeomWasm } from './Dist/ConwayGeomWasm.js'
import { NativeTransform } from '../../src/core/native_types.js'


type NativeVectorGeometryCollection = StdVector<GeometryCollection>
type NativeVectorGeometry = StdVector<GeometryObject>
type NativeVectorMaterial = StdVector<MaterialObject>


export interface StdVector<T> {
  delete(): unknown

  resize(size: number, value?: T): void

  push_back(value: T): void

  size(): number

  empty(): boolean

  set(index: number, value?: T): void
  get(index: number): T

  clear(): void
}


export interface GeometryCollection {

  addComponentWithTransform(geometry: GeometryObject, transform: any): void

  materialIndex: number
  hasDefaultMaterial: boolean
  readonly currentSize: number
}

export interface GeometryObject {
  GetVertexData: () => any
  GetPoint(parameter: number): Vector3
  NormalizeInPlace(): void
  GetVertexDataSize: () => number
  GetIndexData: () => any
  GetIndexDataSize: () => number
  getAllocationSize(): number
  appendGeometry(parameter: GeometryObject): void
  addComponentTransform(transform: any): void
  appendWithTransform(geometry: GeometryObject, transform: any): void
  addComponent(parameter: GeometryObject): void
  clone(): GeometryObject
  applyTransform(parameter: any): void
  min: Vector3
  max: Vector3
  normalized: boolean
  toObj(expressID: number): void
  delete(): void
}


export interface ParamsGetBSplineCurve {
  degree: number
  points2: StdVector< Vector2 >
  points3: StdVector< Vector3 >
  knots: StdVector< number >
  weights: StdVector< number >
}


export interface BSplineSurface {
  active: boolean
  uDegree: number
  vDegree: number
  closedU: boolean
  closedV: boolean
// CurveType has been left out deliberately, it's just metadata -- CS
// weights have been left out, only weightpoints is set from here -- CS
  controlPoints: StdVector< StdVector< Vector3 > >
  uMultiplicity: StdVector< number >
  vMultiplicity: StdVector< number >
  uKnots: StdVector< number >
  vKnots: StdVector< number >
  weightPoints: StdVector< StdVector< number > >
}

export interface CylinderSurface {
  active: boolean
  radius: number
}

export interface IfcProfile3D {
  type: string
  curve: CurveObject
  isConvex: boolean
}

export interface RevolutionSurface {
  active: boolean
  direction: NativeTransform
  profile: IfcProfile3D
}

export interface ExtrusionSurface {
  active: boolean
  direction: Vector3
  profile: ProfileObject
  length: number
}

export interface SurfaceObject {

  transformation: NativeTransform
  bspline: BSplineSurface
  cylinder: CylinderSurface
  revolution: RevolutionSurface
  extrusion: ExtrusionSurface

  normal(): Vector3

  delete(): void
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
  delete(): void
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

export interface ParamsAxis1Placement3D {
  position: any
  zAxisRef: any
  normalizeZ: boolean
}

export interface ParamsGetBooleanResult {
  flatFirstMesh: any
  flatSecondMesh: any
  operatorType: number
}

export interface ParamsRelVoidSubtract {
  flatFirstMesh: any
  flatSecondMesh: any
  operatorType: number
  parentMatrix: any
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
   * @param initialSize number - initial size of the vector (optional)
   * @return {NativeVectorGeometry} - a native std::vector<GeometryObject> from the wasm module
   */
  nativeVectorGeometry(initialSize?: number): NativeVectorGeometry {
    const nativeVectorGeometry_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.geometryArray)()) as NativeVectorGeometry

    if (initialSize) {
      const defaultGeometry = (new (this.wasmModule.IfcGeometry)) as GeometryObject
      // resize has a required second parameter to set default values
      nativeVectorGeometry_.resize(initialSize, defaultGeometry)
    }

    return nativeVectorGeometry_
  }

  /**
   *
   * @return {NativeVectorGeometry} - a native std::vector<GeometryObject> from the wasm module
   */
  nativeVectorVectorDouble(): StdVector< StdVector< number > > {
    const nativeVectorVectorDouble_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.vectorVectorDouble)()) as StdVector< StdVector< number > >

    return nativeVectorVectorDouble_
  }

  /**
   *
   * @param initialSize number - initial size of the vector (optional)
   * @return {NativeVectorGeometry} - a native std::vector<GeometryObject> from the wasm module
   */
  nativeVectorDouble(initialSize?: number): StdVector< number > {
    const nativeVectorDouble_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.vectorDouble)()) as StdVector< number >

    if (initialSize !== void 0) {
      // resize has a required second parameter to set default values
      nativeVectorDouble_.resize(initialSize, 0)
    }

    return nativeVectorDouble_
  }

  /**
   * Create a native geometry collection.
   *
   * @return {GeometryCollection}
   */
  nativeGeometryCollection(): GeometryCollection {
    const nativeGeometryCollection =
      (new (this.wasmModule.IfcGeometryCollection)()) as GeometryCollection

    return nativeGeometryCollection
  }

  /**
   *
   * @return {NativeVectorGeometry} - a native std::vector<GeometryObject> from the wasm module
   */
  nativeVectorGeometryCollection(): NativeVectorGeometryCollection {
    const nativeVectorGeometryCollection_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.geometryCollectionArray)()) as NativeVectorGeometryCollection

    return nativeVectorGeometryCollection_
  }

  /**
   *
   * @param from The material to create the native material from
   * @return {MaterialObject} The created canonical material.
   */
  nativeMaterial(from: CanonicalMaterial): MaterialObject {
    const native: MaterialObject = {

      alphaCutoff: 0,
      alphaMode: toAlphaMode(this.wasmModule, from.blend),
      base: {
        x: from.baseColor[0],
        y: from.baseColor[1],
        z: from.baseColor[2],
        w: from.baseColor[3],
      },
      doubleSided: from.doubleSided,
      /* eslint-disable no-magic-numbers */
      ior: from.ior ?? 1.4,
      metallic: from.metalness ?? 1.0,
      roughness: from.roughness ?? 1.0,
      specular: from.specular !== void 0 ? {
        x: from.specular[0],
        y: from.specular[1],
        z: from.specular[2],
        w: from.specular[3],
      } : void 0,
    }
    /* eslint-enable no-magic-numbers */
    return native
  }

  /**
   *
   * @param initialSize number - initial size of the vector (optional)
   * @return {NativeVectorMaterial} - a native std::vector<MaterialObject> from the wasm module
   */
  nativeVectorMaterial(initialSize?: number): NativeVectorMaterial {
    const nativeVectorMaterial_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.materialArray)()) as NativeVectorMaterial

    if (initialSize) {
      // resize has a required second parameter to set default values
      nativeVectorMaterial_.resize(initialSize)
    }

    return nativeVectorMaterial_
  }

  /**
   *
   * @return {Promise<boolean>} - initialization status
   */
  async initialize(fileHandler?: LocateFileHandlerFn): Promise<boolean> {
    if (this.wasmModule === void 0) {
      // eslint-disable-next-line new-cap
      this.wasmModule = await ConwayGeomWasm({ noInitialRun: true, locateFile: fileHandler })
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
   * Get a B-Spline Curve
   *
   * @param parameters
   * @return {CurveObject} - The native curve object.
   */
  getBSplineCurve(parameters: ParamsGetBSplineCurve ) : CurveObject {
    const result = this.wasmModule.getBSplineCurve( parameters )
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
   * @param parameters
   * @return {GeometryObject}
   */
  relVoidSubtract(parameters: ParamsRelVoidSubtract): GeometryObject {
    const result = this.wasmModule.relVoidSubtract(parameters)
    return result
  }

  /**
   * Convert geometry to gltf.
   *
   * @param mat1
   * @param mat2
   * @return {any} matrix result of the multiplication
   */
  multiplyNativeMatrices(mat1: NativeTransform, mat2: NativeTransform): any {
    const result = this.wasmModule.multiplyNativeMatrices(mat1, mat2)
    return result
  }

  /**
   *
   * @param geometry Vector of native geometry collection objects
   * @param materials Vector of native materials indexed by geometry
   * @param isGlb  boolean if the output should be a single GLB file
   * @param outputDraco boolean should the output use Draco compression
   * @param fileUri string of base name for output files
   * @param geometryOffset The offset into the geometry vector to use to start
   *
   * @return {ResultsGltf} boolean success + buffers + file uris
   */
  toGltf(
      geometry: StdVector<GeometryCollection>,
      materials: StdVector<MaterialObject>,
      isGlb: boolean,
      outputDraco: boolean,
      fileUri: string,
      geometryOffset: number = 0,
      geometryCount: number = geometry.size() ):
      ResultsGltf {
    return this.wasmModule.geometryToGltf(
        geometry,
        materials,
        isGlb,
        outputDraco,
        fileUri,
        geometryOffset,
        geometryCount)
  }

  /**
   *
   * @param parameters - ParamsGetAxis2Placement2D structure
   * @return {any} - native Axis2Placement2D structure
   */
  getAxis1Placement3D(parameters: ParamsAxis1Placement3D): any {
    return this.wasmModule.getAxis1Placement(parameters)
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
   * @return {any} identity matrix
   */
  getIdentityTransform():any {
    return this.wasmModule.getIdentityTransform()
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
