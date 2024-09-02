import { NativeTransform4x4 } from '../native_transform'


/** Parameter set for a local placement */
export interface ParamsLocalPlacement {
  useRelPlacement: boolean
  axis2Placement: NativeTransform4x4
  relPlacement: NativeTransform4x4
}
