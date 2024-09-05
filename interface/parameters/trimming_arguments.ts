import { TrimmingSelect } from './trimming_select'


/** Parameters for trimming curves */
export interface TrimmingArguments {
  exist: boolean
  start: TrimmingSelect | undefined
  end: TrimmingSelect | undefined
}
