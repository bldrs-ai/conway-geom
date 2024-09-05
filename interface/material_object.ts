import { ConwayGeometryWasm } from './conway_geometry'
import { Vector4 } from './vector4'


/**
 * A native code material.
 */
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


/* eslint-disable no-shadow,no-unused-vars,no-magic-numbers */
export enum BlendMode {
  OPAQUE = 0,
  BLEND = 1,
  MASK = 2,
}

/**
 * @return {any}
 */
export function toAlphaMode(wasmModule: ConwayGeometryWasm, blendMode: BlendMode): any {

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
