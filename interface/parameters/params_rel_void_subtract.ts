import { GeometryObject } from '../geometry_object'
import { NativeTransform4x4 } from '../native_transform'
import { StdVector } from '../std_vector'


/** Parameter set to perform a rel void subtract */
export interface ParamsRelVoidSubtract {
  flatFirstMesh: StdVector< GeometryObject >
  flatSecondMesh: StdVector< GeometryObject >
  operatorType: number
  parentMatrix: NativeTransform4x4
}
