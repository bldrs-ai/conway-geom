
import { Bound3DObject } from './bound_3D_object'
import { CurveObject } from './curve_object'
import { GeometryCollection } from './geometry_collection'
import { GeometryObject } from './geometry_object'
import { MaterialObject } from './material_object'
import { ParamsGetLoop } from './parameters/params_get_loop'
import { StdVector } from './std_vector'
import { ParamsAddFaceToGeometry } from './parameters/params_add_face_to_geometry'
import {
  ParamsCartesianTransformationOperator3D,
} from './parameters/params_cartesian_transform_operator_3D'
import { ProfileObject } from './profile_object'
import { ParamsTransformProfile } from './parameters/params_transform_profile'
import { ParamsGetPolyCurve } from './parameters/params_get_poly_curve'
import { ParamsPolygonalFaceSet } from './parameters/params_polygonal_face_set'
import {
  ParamsGetTriangulatedFaceSetGeometry,
} from './parameters/params_get_triangulated_face_set_geometry'
import { ParamsGetCShapeCurve } from './parameters/params_get_c_shape_curve'
import { ParamsGetIShapeCurve } from './parameters/params_get_i_shape_curve'
import { ParamsGetLShapeCurve } from './parameters/params_get_l_shape_curve'
import { ParamsGetTShapeCurve } from './parameters/params_get_t_shape_curve'
import { ParamsGetUShapeCurve } from './parameters/params_get_u_shape_curve'
import { ParamsGetZShapeCurve } from './parameters/params_get_z_shape_curve'
import { ParamsGetIfcCircle } from './parameters/params_get_ifc_circle'
import { ParamsGetIfcLine } from './parameters/params_get_ifc_line'
import { ParamsGetBSplineCurve } from './parameters/params_get_bspline_curve'
import { ParamsGetIfcIndexedPolyCurve } from './parameters/params_get_ifc_indexed_poly_curve'
import { ParamsGetCircleCurve } from './parameters/params_get_circle_curve'
import { ParamsGetEllipseCurve } from './parameters/params_get_ellipse_curve'
import { ParamsCreateNativeIfcProfile } from './parameters/params_create_native_ifc_profile'
import { ParamsCreateBound3D } from './parameters/params_create_bound_3D'
import { ParamsGetHalfspaceSolid } from './parameters/params_get_half_space_solid'
import { ParamsGetRectangleProfileCurve } from './parameters/params_get_rectangle_profile_curve'
import {
  ParamsGetPolygonalBoundedHalfspace,
} from './parameters/params_get_polygonal_bounded_halfspace'
import { ParamsGetExtrudedAreaSolid } from './parameters/params_get_extruded_area_solid'
import { ParamsGetBooleanResult } from './parameters/params_get_boolean_result'
import { ParamsRelVoidSubtract } from './parameters/params_rel_void_subtract'
import { NativeTransform3x3, NativeTransform4x4 } from './native_transform'
import { ResultsGltf } from './results_gltf'
import { ParamsAxis1Placement3D } from './parameters/params_axis_1_placement_3D'
import { ParamsGetAxis2Placement2D } from './parameters/params_get_axis_2_placement_2D'
import { ParamsAxis2Placement3D } from './parameters/params_axis_2_placement_3D'
import { ParamsLocalPlacement } from './parameters/params_local_placement'
import { ParamsGetSweptDiskSolid } from './parameters/params_get_swept_disk_solid'
import { ParseBuffer } from './parse_buffer'


let ConwayGeomWasm: any

/**
 * Load the WebAssembly module based on the environment
 */
async function loadWasmModule() {
  if (process.env.PLATFORM === 'web') {
    // Load browser-specific WebAssembly module
    const module = await import('../Dist/ConwayGeomWasmWeb.js')
    ConwayGeomWasm = module.default
  } else {
    // Load Node.js-specific WebAssembly module
    const module = await import('../Dist/ConwayGeomWasmNode.js')
    ConwayGeomWasm = module.default
  }
}

export default ConwayGeomWasm


export type ConwayGeometryWasm = typeof ConwayGeomWasm

export type FileHandlerFunction = (path: string, prefix: string) => string

/**
 * Internal interface for wasm module, geometry processing
 * OBJ + GLTF + GLB (Draco) Conversions
 */
export class ConwayGeometry {
  public wasmModule?: ConwayGeometryWasm
  initialized = false

  private parseBuffers_: ParseBuffer[] = []

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
   * @return {StdVector< GeometryObject >} - a native std::vector<GeometryObject> from the
   * wasm module
   */
  nativeVectorGeometry(initialSize?: number): StdVector< GeometryObject > {
    const nativeVectorGeometry_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.geometryArray)()) as StdVector< GeometryObject >

    if (initialSize) {
      const defaultGeometry = (new (this.wasmModule.IfcGeometry)) as GeometryObject
      // resize has a required second parameter to set default values
      nativeVectorGeometry_.resize(initialSize, defaultGeometry)
    }

    return nativeVectorGeometry_
  }

  /**
   *
   * @return {GeometryObject} - an empty native geometry object
   */
  nativeGeometry(initialSize?: number): GeometryObject {

    const nativeGeometry = (new (this.wasmModule.IfcGeometry)) as GeometryObject

    return nativeGeometry
  }

  /**
   *
   *
   * @return {ParseBuffer} - a parse buffer
   */
  nativeParseBuffer(): ParseBuffer {

    if ( this.parseBuffers_.length > 0 ) {

      return this.parseBuffers_.pop()!
    }

    const nativeParseBuffer =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.ParseBuffer)()) as ParseBuffer

    return nativeParseBuffer
  }

  /**
   *
   *
   * @param buffer The parse buffer to free
   */
  freeParseBuffer( buffer: ParseBuffer ): void {

    this.parseBuffers_.push( buffer )
  }

  /**
   * Allocate a native vector of vectors of doubles.
   *
   * @return {StdVector< GeometryObject >} - a native std::vector<GeometryObject> from the
   * wasm module
   */
  nativeVectorVectorDouble(): StdVector<StdVector<number>> {
    const nativeVectorVectorDouble_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.vectorVectorDouble)()) as StdVector<StdVector<number>>

    return nativeVectorVectorDouble_
  }

  /**
   * Allocate a native vector of doubles.
   *
   * @param initialSize number - initial size of the vector (optional)
   * @return {StdVector< GeometryObject >} - a native std::vector<GeometryObject> from the
   * wasm module
   */
  nativeVectorDouble(initialSize?: number): StdVector<number> {
    const nativeVectorDouble_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.vectorDouble)()) as StdVector<number>

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
   * @return {StdVector< GeometryObject >} - a native std::vector<GeometryObject> from the
   * wasm module
   */
  nativeVectorGeometryCollection(): StdVector< GeometryCollection > {
    const nativeVectorGeometryCollection_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.geometryCollectionArray)()) as StdVector< GeometryCollection >

    return nativeVectorGeometryCollection_
  }

  /**
   *
   * @param initialSize number - initial size of the vector (optional)
   * @return {StdVector< MaterialObject >} - a native std::vector<MaterialObject> from the
   * wasm module
   */
  nativeVectorMaterial(initialSize?: number): StdVector< MaterialObject > {
    const nativeVectorMaterial_ =
      // eslint-disable-next-line new-cap
      (new (this.wasmModule.materialArray)()) as StdVector< MaterialObject >

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
  async initialize(fileHandler?: FileHandlerFunction): Promise<boolean> {
    // Wait for the WebAssembly module to load if it's not already set
    if (!ConwayGeomWasm) {
      await loadWasmModule();
    }

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
   * @param parameters - ParamsTransformProfile parsed from data model
   * @return {ProfileObject} - Native Profile object
   */
  transformProfile(parameters: ParamsTransformProfile): ProfileObject {
    const result = this.wasmModule.transformProfile(parameters)
    return result
  }

  /**
   *
   * @param parameters - ParamsGetPolyCurve parsed from data model
   * @return {CurveObject} - Native Curve object
   */
  getPolyCurve(parameters: ParamsGetPolyCurve):CurveObject {
    const result = this.wasmModule.getPolyCurve(parameters)
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
   * @param parameters ParamsGetPolygonalFaceSetGeometry parsed from data model
   * @return {GeometryObject} - Native Geometry Object
   */
  getTriangulatedFaceSetGeometry(parameters:ParamsGetTriangulatedFaceSetGeometry):GeometryObject {
    const result = this.wasmModule.getTriangulatedFaceSetGeometry(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetCShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getCShapeCurve(parameters: ParamsGetCShapeCurve): CurveObject {
    const result = this.wasmModule.getCShapeCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetIShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getIShapeCurve(parameters: ParamsGetIShapeCurve): CurveObject {
    const result = this.wasmModule.getIShapeCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetLShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getLShapeCurve(parameters: ParamsGetLShapeCurve): CurveObject {
    const result = this.wasmModule.getLShapeCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetTShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getTShapeCurve(parameters: ParamsGetTShapeCurve): CurveObject {
    const result = this.wasmModule.getTShapeCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetUShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getUShapeCurve(parameters: ParamsGetUShapeCurve): CurveObject {
    const result = this.wasmModule.getUShapeCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetZShapeCurve parsed from data model
   * @return {CurveObject}
   */
  getZShapeCurve(parameters: ParamsGetZShapeCurve): CurveObject {
    const result = this.wasmModule.getZShapeCurve(parameters)
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
   * @param parameters
   * @return {CurveObject} - Native Curve Object
   */
  getIfcLine(parameters: ParamsGetIfcLine): CurveObject {
    const result = this.wasmModule.getIfcLine(parameters)
    return result
  }

  /**
   * Get a B-Spline Curve
   *
   * @param parameters
   * @return {CurveObject} - The native curve object.
   */
  getBSplineCurve(parameters: ParamsGetBSplineCurve): CurveObject {
    const result = this.wasmModule.getBSplineCurve(parameters)
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
   * @param parameters ParamsGetCirlceCurve parsed from data model
   * @return {CurveObject}
   */
  getCircleHoleCurve(parameters: ParamsGetCircleCurve): CurveObject {
    const result = this.wasmModule.getCircleHoleCurve(parameters)
    return result
  }

  /**
   *
   * @param parameters ParamsGetEllipseCurve parsed from data model
   * @return {CurveObject}
   */
  getEllipseCurve(parameters: ParamsGetEllipseCurve): CurveObject {
    const result = this.wasmModule.getEllipseCurve(parameters)
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
   * @param parameters 
   * @returns {GeometryObject}
   */
  getSweptDiskSolid(parameters: ParamsGetSweptDiskSolid): GeometryObject {
    const result = this.wasmModule.getSweptDiskSolid(parameters);
    return result;
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
   * @param parameters ParamsGetRectangleProfileCurve parsed from data model
   * @return {CurveObject}
   */
  getRectangleHollowProfileHole(parameters: ParamsGetRectangleProfileCurve): CurveObject {
    const result = this.wasmModule.getRectangleHollowProfileHole(parameters)
    return result
  }

  /**
   *
   * @param parameters
   * @return {GeometryObject}
   */
  getPolygonalBoundedHalfspace(parameters:ParamsGetPolygonalBoundedHalfspace):GeometryObject {
    const result = this.wasmModule.getPolygonalBoundedHalfspace(parameters)
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
  multiplyNativeMatrices(mat1: NativeTransform4x4, mat2: NativeTransform4x4): NativeTransform4x4 {
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
   * @return {ResultsGltf} boolean success + buffers + file uris
   */
  toGltf(
      geometry: StdVector<GeometryCollection>,
      materials: StdVector<MaterialObject>,
      isGlb: boolean,
      outputDraco: boolean,
      fileUri: string,
      geometryOffset: number = 0,
      geometryCount: number = geometry.size()):
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
  getAxis1Placement3D(parameters: ParamsAxis1Placement3D): NativeTransform4x4 {
    return this.wasmModule.getAxis1Placement(parameters)
  }

  /**
   *
   * @param parameters - ParamsGetAxis2Placement2D structure
   * @return {any} - native Axis2Placement2D structure
   */
  getAxis2Placement2D(parameters: ParamsGetAxis2Placement2D): NativeTransform3x3 {
    return this.wasmModule.getAxis2Placement2D(parameters)
  }

  /**
   *
   * @param parameters - ParamsAxis2Placement3D structure
   * @return {any} - native Axis2Placement3D structure
   */
  getAxis2Placement3D(parameters: ParamsAxis2Placement3D): NativeTransform4x4 {
    return this.wasmModule.getAxis2Placement3D(parameters)
  }

  /**
   *
   * @param parameters - ParamsLocalPlacement structure
   * @return {any} = native LocalPlacement structure
   */
  getLocalPlacement(parameters: ParamsLocalPlacement): NativeTransform4x4 {
    return this.wasmModule.getLocalPlacement(parameters)
  }

  /**
   *
   * @return {any} identity matrix
   */
  getIdentityTransform(): NativeTransform4x4 {
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
