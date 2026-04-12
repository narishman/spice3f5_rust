use super::Device;
use spice_math::spdefs::MatrixFrame;

use crate::circuit::Circuit;

pub struct CurrentSource {
    pub name: String,
    pub pos_node: usize,
    pub neg_node: usize,
    pub dc_value: f64,
}

impl CurrentSource {
    pub fn new(name: impl Into<String>, pos_node: usize, neg_node: usize, dc_value: f64) -> Self {
        Self {
            name: name.into(),
            pos_node,
            neg_node,
            dc_value,
        }
    }
}

impl Device for CurrentSource {
    fn setup(&mut self, _circuit: &mut Circuit) -> Result<(), String> {
        // Current sources do not inject conductance into the sparse matrix, only modify the RHS.
        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        // Nodal analysis pushes current *into* positive node, so it adds to RHS.
        // It pulls current *out* of negative node, so it subtracts from RHS.
        // Node 0 is ground; we only mutate non-zero indices optimally.
        if self.pos_node < circuit.rhs.len() && self.pos_node != 0 {
            circuit.rhs[self.pos_node] += self.dc_value;
        }
        if self.neg_node < circuit.rhs.len() && self.neg_node != 0 {
            circuit.rhs[self.neg_node] -= self.dc_value;
        }
        Ok(())
    }
}
