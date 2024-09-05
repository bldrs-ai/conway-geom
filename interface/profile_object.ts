import { CurveObject } from './curve_object'
import { Deletable } from './deletable'
import { StdVector } from './std_vector'
import { WavefrontDumpable } from './wavefront_dumpable'


/**
 * A profile object (representing ap profile outline, which includes a main outer
 * curve and inner hole curves, as well as potentially multiple child profiles.
 */
export interface ProfileObject extends Deletable, WavefrontDumpable {
  getType: () => string
  getCurve: () => CurveObject
  getHoles: () => StdVector< CurveObject >
  isConvex: () => boolean
  isComposite: () => boolean
  getProfiles: () => StdVector< ProfileObject >
}
