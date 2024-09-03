import { Deletable } from '../deletable'
import { NativeTransform } from '../native_transform'
import { ProfileObject } from '../profile_object'


/** Parameter set for transforming a profile */
export interface ParamsTransformProfile extends Deletable {
  transformation: NativeTransform // glm::dmat3
  profile: ProfileObject
}
