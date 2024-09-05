import { StdVector } from './std_vector'


/** Results from gltf  */
export interface ResultsGltfBase {
  success: boolean
  bufferUris?: StdVector< string >
  buffers?: StdVector< StdVector< number > >
}

export interface ResultsGltfSuccess extends ResultsGltfBase {
  success: true
  bufferUris: StdVector< string >
  buffers: StdVector< StdVector< number > >
}

export interface ResultsGltfFailure extends ResultsGltfBase {
  success: false
  bufferUris: undefined
  buffers: undefined
}

export type ResultsGltf = ResultsGltfSuccess | ResultsGltfFailure
