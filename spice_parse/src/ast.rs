#[derive(Debug, PartialEq, Clone)]
pub enum Stmt {
    Device(DeviceStmt),
    Control(ControlStmt),
    Model(ModelStmt),
    Comment(String),
}

#[derive(Debug, PartialEq, Clone)]
pub struct DeviceStmt {
    pub ident: String,              // e.g., "R1"
    pub nodes: Vec<String>,         // e.g., ["in", "out"]
    pub model: Option<String>,      // Optional model name for D, Q, M devices
    pub value: f64,                 // e.g., 10000.0 from "10k"
    pub params: Vec<(String, f64)>, // Params like `w=10u` or trailing numeric pairs.
}

#[derive(Debug, PartialEq, Clone)]
pub enum ControlStmt {
    Op,
    Dc {
        source: String,
        start: f64,
        stop: f64,
        step: f64,
    },
    Tran {
        step: f64,
        stop: f64,
    },
    End,
}

#[derive(Debug, PartialEq, Clone)]
pub struct ModelStmt {
    pub name: String,
    pub model_type: String,         // e.g. "NPN", "PNP", "NMOS", "PMOS", "D"
    pub params: Vec<(String, f64)>, // e.g. BF=100 IS=1e-14
}
