import { Deletable } from '../deletable'
import { GeometryObject } from '../geometry_object'
import { StdVector } from '../std_vector'


/** Parameters to get boolean result */
export interface ParamsGetBooleanResult extends Deletable {
  flatFirstMesh: StdVector< GeometryObject >
  flatSecondMesh: StdVector< GeometryObject >
  operatorType: number
}
