import { NativeTransform } from '../native_transform'
import { ProfileObject } from '../profile_object'


/** Parameter object for a revolution surface */
export interface RevolutionSurface {
  active: boolean
  direction: NativeTransform
  profile: ProfileObject
}
