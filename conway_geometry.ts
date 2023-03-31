let ConwayWasm: any;

//@ts-ignore
if (typeof self !== 'undefined' && self.crossOriginIsolated) {
    ConwayWasm = require("./conway_geom_wasm-mt");
}
else {
    ConwayWasm = require("./conway_geom_wasm");
}


export interface IfcGeometry {
    GetVertexData(): number;
    GetVertexDataSize(): number;
    GetIndexData(): number;
    GetIndexDataSize(): number;
}

export class ConwayGeometry {

}

