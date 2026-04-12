#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeType {
    Voltage,
    Current,
}

#[derive(Debug, Clone)]
pub struct Node {
    pub name: String,
    pub node_type: NodeType,
    pub number: usize, // Internal matrix index
    pub ic: f64,
    pub nodeset: f64,
    pub ic_given: bool,
    pub ns_given: bool,
}

impl Node {
    pub fn new(name: impl Into<String>, node_type: NodeType, number: usize) -> Self {
        Self {
            name: name.into(),
            node_type,
            number,
            ic: 0.0,
            nodeset: 0.0,
            ic_given: false,
            ns_given: false,
        }
    }
}
