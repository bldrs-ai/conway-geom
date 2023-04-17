
import { default as conway_geom_wasm } from "./conway_geom_wasm.js"



export interface GeometryObject {
    GetVertexData: () => any
    GetVertexDataSize: () => number
    GetIndexData: () => any
    GetIndexDataSize: () => number
    AddGeometry(parameter: GeometryObject): void
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

export type LocateFileHandlerFn = (path: string, prefix: string) => string

export class ConwayGeometry {
    modelId: number = -1
    //private module: any

    public wasmModule: undefined | any = undefined
    fs: undefined | any = undefined
    wasmPath: string = ""
    isWasmPathAbsolute = false
    initialized = false

    constructor() {
    }

    async initialize() {
        if (!this.wasmModule) {
            this.wasmModule = await new conway_geom_wasm()
        }
        //console.log(this.wasmModule)
        this.modelId = this.wasmModule.InitializeGeometryProcessor()

        //console.log("modelId: " + this.modelId)

        this.initialized = true
    }

    addGeometry(parameter: GeometryObject) {
        this.wasmModule.AddGeometry(parameter)
    }

    getGeometry(parameters: ParamsPolygonalFaceSet): GeometryObject {
        const result = this.wasmModule.GetGeometry(this.modelId, parameters)
        return result
    }

    toGltf(geometry: GeometryObject, isGlb: boolean, outputDraco: boolean, fileUri: string): ResultsGltf {
        return this.wasmModule.GeometryToGltf(this.modelId, geometry, isGlb, outputDraco, fileUri)
    }

    toObj(geometry: GeometryObject): string {
        return this.wasmModule.GeometryToObj(this.modelId, geometry, 0)
    }

    destroy() {
        this.wasmModule.FreeGeometryProcessor(this.modelId)
        this.modelId = -1
        this.initialized = false;
    }
}

