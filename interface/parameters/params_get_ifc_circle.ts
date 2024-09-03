import { NativeTransform3x3, NativeTransform4x4 } from '../native_transform'
import { ParamsGetIfcTrimmedCurve } from './params_get_ifc_trimmed_curve'


/** Parameter set to get an ifc circle curve */
export interface ParamsGetIfcCircle {
  dimensions: number
  axis2Placement2D: NativeTransform3x3
  axis2Placement3D: NativeTransform4x4
  radius: number
  radius2: number
  paramsGetIfcTrimmedCurve: ParamsGetIfcTrimmedCurve
}
