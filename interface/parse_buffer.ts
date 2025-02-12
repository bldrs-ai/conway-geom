import { Deletable } from './deletable'

/**
 * Wrapper interface for a native std::vector as exposed by embind.
 */
export interface ParseBuffer extends Deletable {

  resize(size: number): number

  size(): number

  capacity(): number

  data(): number
}
