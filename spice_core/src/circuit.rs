use crate::node::{Node, NodeType};
use spice_math::spdefs::MatrixFrame;
use std::collections::HashMap;

pub struct Circuit {
    pub matrix: MatrixFrame,
    pub nodes: Vec<Node>,
    pub node_map: HashMap<String, usize>, // Maps name to node index

    // States
    // In SPICE, there are 8 state vectors used for integration history.
    pub states: Vec<Vec<f64>>,
    pub rhs_old: Vec<f64>,
    pub rhs: Vec<f64>,

    // Iteration & Tolerances
    pub reltol: f64,
    pub abstol: f64,
    pub volt_tol: f64,
    pub dc_max_iter: usize,
    pub tran_max_iter: usize,
    pub gmin: f64,
    pub bypass: bool,

    // Convergence and integration status
    pub ni_state: i32,
    pub noncon: usize, // Number of non-convergent nodes

    // Transient & Active States
    pub num_states: usize,   // Number of active variable states allocated
    pub ag: [f64; 7],        // Integration coefficients
    pub delta: f64,          // Current timestep
    pub delta_old: [f64; 7], // History of timesteps
}

impl Circuit {
    pub fn new() -> Self {
        let mut d = Self {
            matrix: MatrixFrame::new(0, false), // Start with matrix dimension 0
            nodes: Vec::new(),
            node_map: HashMap::new(),
            states: vec![Vec::new(); 8],
            rhs_old: Vec::new(),
            rhs: Vec::new(),
            reltol: 1e-3,
            abstol: 1e-12,
            volt_tol: 1e-6,
            dc_max_iter: 100,
            tran_max_iter: 10,
            gmin: 1e-12,
            bypass: false,
            ni_state: 0,
            noncon: 0,
            num_states: 0,
            ag: [0.0; 7],
            delta: 0.0,
            delta_old: [0.0; 7],
        };

        // Create the common Ground node (0)
        d.nodes.push(Node::new("0", NodeType::Voltage, 0));
        d.node_map.insert("0".to_string(), 0);

        d
    }

    /// Add or return an existing node by name
    pub fn get_or_create_node(&mut self, name: &str, node_type: NodeType) -> usize {
        if let Some(&idx) = self.node_map.get(name) {
            return idx;
        }

        let new_idx = self.nodes.len();
        self.nodes.push(Node::new(name, node_type, new_idx));
        self.node_map.insert(name.to_string(), new_idx);

        // Expand vectors accommodating new nodes
        self.rhs.resize(new_idx + 1, 0.0);
        self.rhs_old.resize(new_idx + 1, 0.0);

        new_idx
    }
}

impl Default for Circuit {
    fn default() -> Self {
        Self::new()
    }
}
