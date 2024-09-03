import { Deletable } from './deletable'

/**
 * Wrapper interface for a native std::vector as exposed by embind.
 */
export interface StdVector<T> extends Deletable {

  resize(size: number, value?: T): void

  push_back(value: T): void

  size(): number

  empty(): boolean

  set(index: number, value?: T): void
  get(index: number): T

  clear(): void
}
