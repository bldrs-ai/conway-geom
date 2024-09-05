/** Interface for dumping a wavefront OBJ file from an object */
export interface WavefrontDumpable {

  dumpToOBJ: ( preamble: string ) => string
}
