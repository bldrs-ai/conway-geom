import { CurveObject } from '../curve_object'


/** Parameters for creating a 3D bound */
export interface ParamsCreateBound3D {
  curve: CurveObject // conway::geometry::IfcCurve
  orientation: boolean
  type: number // uint32_t

  /**
   * This loop retraces itself - every edge traversed exactly twice, once in
   * each direction - so it encloses no area and cannot trim. Per ISO
   * 10303-42 that means the face covers the whole surface, the same as a
   * VERTEX_LOOP. Decided by the front end from the ORIENTED_EDGEs; see
   * bldrs-ai/conway#595.
   */
  seam: boolean

  /**
   * This loop CONTAINS a retraced edge pair - one edge walked once in each
   * direction - so it wraps the seam of a surface that is closed in one
   * parameter, and the face covers that parameter's full period. Distinct
   * from `seam` above, which is the whole loop retracing. Decided by the
   * front end from the ORIENTED_EDGEs; see bldrs-ai/conway#611.
   */
  seamPair: boolean
}
