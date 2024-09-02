import { Deletable } from './deletable'

/**
 * A native transform (read matrix).
 *
 * Note, is expected from context to know what size is required for "setvalues"
 * as the matrix might be 2x2, 3x3 or 4x4. If this needs to be determined at runtime,
 * getValues with a .length can be used to determine the number of elements and
 * corresponding matrix.
 */
export interface NativeTransform extends Deletable {

  /**
   * Get the values for this transform matrix.
   *
   * @return {Readonly<number[]>} The values in this matrix.
   */
  getValues(): Readonly<number[]>

  /**
   *  Set the values in this matrix.
   *
   * @param values The values to set.
   */
  setValues( values: ArrayLike<number> ): void
}

export interface NativeTransform3x3 extends NativeTransform {

}

export interface NativeTransform4x4 extends NativeTransform {

}
