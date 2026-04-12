use super::Device;
use crate::circuit::Circuit;
use crate::node::NodeType;
use spice_math::spdefs::MatrixFrame;

pub struct VoltageSource {
    pub name: String,
    pub pos_node: usize,
    pub neg_node: usize,
    pub dc_value: f64,
    pub branch_node: usize,

    pub pos_ibr_ptr: Option<usize>,
    pub neg_ibr_ptr: Option<usize>,
    pub ibr_pos_ptr: Option<usize>,
    pub ibr_neg_ptr: Option<usize>,
}

impl VoltageSource {
    pub fn new(name: impl Into<String>, pos_node: usize, neg_node: usize, dc_value: f64) -> Self {
        Self {
            name: name.into(),
            pos_node,
            neg_node,
            dc_value,
            branch_node: 0,
            pos_ibr_ptr: None,
            neg_ibr_ptr: None,
            ibr_pos_ptr: None,
            ibr_neg_ptr: None,
        }
    }
}

impl Device for VoltageSource {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        if self.branch_node == 0 {
            // MNA Voltage Sources require a dedicated current branch equation.
            let branch_name = format!("branch_{}", self.name);
            self.branch_node = circuit.get_or_create_node(&branch_name, NodeType::Current);
        }

        self.pos_ibr_ptr = circuit.matrix.get_element(self.pos_node, self.branch_node);
        self.neg_ibr_ptr = circuit.matrix.get_element(self.neg_node, self.branch_node);
        self.ibr_pos_ptr = circuit.matrix.get_element(self.branch_node, self.pos_node);
        self.ibr_neg_ptr = circuit.matrix.get_element(self.branch_node, self.neg_node);
        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        // Voltage Source adds elements representing I(Vsrc) logic.
        if let Some(idx) = self.pos_ibr_ptr {
            circuit.matrix.elements_arena[idx].real += 1.0;
        }
        if let Some(idx) = self.neg_ibr_ptr {
            circuit.matrix.elements_arena[idx].real -= 1.0;
        }
        if let Some(idx) = self.ibr_pos_ptr {
            circuit.matrix.elements_arena[idx].real += 1.0;
        }
        if let Some(idx) = self.ibr_neg_ptr {
            circuit.matrix.elements_arena[idx].real -= 1.0;
        }

        // Setup the right hand side with the static DC voltage value
        if self.branch_node < circuit.rhs.len() {
            circuit.rhs[self.branch_node] += self.dc_value;
        }
        Ok(())
    }
}
