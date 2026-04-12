use crate::circuit::Circuit;
use crate::devices::Device;
use crate::integrate::integrate;
use crate::node::NodeType;

pub struct Inductor {
    pub name: String,
    pub pos_node: usize,
    pub neg_node: usize,
    pub inductance: f64,
    pub branch_node: usize,
    pub state_idx: usize, // Pointer to state buffer slot

    pub ibr_ibr_ptr: Option<usize>,
    pub pos_ibr_ptr: Option<usize>,
    pub neg_ibr_ptr: Option<usize>,
    pub ibr_pos_ptr: Option<usize>,
    pub ibr_neg_ptr: Option<usize>,
}

impl Inductor {
    pub fn new(name: impl Into<String>, pos_node: usize, neg_node: usize, inductance: f64) -> Self {
        Self {
            name: name.into(),
            pos_node,
            neg_node,
            inductance,
            branch_node: 0,
            state_idx: 0,
            ibr_ibr_ptr: None,
            pos_ibr_ptr: None,
            neg_ibr_ptr: None,
            ibr_pos_ptr: None,
            ibr_neg_ptr: None,
        }
    }
}

impl Device for Inductor {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        self.state_idx = circuit.num_states;
        circuit.num_states += 2; // Flux (Phi) and Voltage (V)

        if self.branch_node == 0 {
            let branch_name = format!("branch_{}", self.name);
            self.branch_node = circuit.get_or_create_node(&branch_name, NodeType::Current);
        }

        self.ibr_ibr_ptr = circuit
            .matrix
            .get_element(self.branch_node, self.branch_node);
        self.pos_ibr_ptr = circuit.matrix.get_element(self.pos_node, self.branch_node);
        self.neg_ibr_ptr = circuit.matrix.get_element(self.neg_node, self.branch_node);
        self.ibr_pos_ptr = circuit.matrix.get_element(self.branch_node, self.pos_node);
        self.ibr_neg_ptr = circuit.matrix.get_element(self.branch_node, self.neg_node);

        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        if circuit.delta > 0.0 {
            // Transient integration: compute flux
            let current = circuit.rhs_old[self.branch_node]; // I_L
            circuit.states[0][self.state_idx] = self.inductance * current; // Flux Phi = L * I

            let mut req = 0.0;
            let mut veq = 0.0;

            // integrate computes: req = L * 2/h, veq = -(2/h * Phi_1 + V_1)
            integrate(circuit, &mut req, &mut veq, self.inductance, self.state_idx);

            // MNA formulation for Inductor branch: V_pos - V_neg - L * di/dt = 0
            // Since L*di/dt = req * i + veq
            // V_pos - V_neg - req * i = veq

            if let Some(idx) = self.ibr_ibr_ptr {
                circuit.matrix.elements_arena[idx].real -= req;
            }
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

            circuit.rhs[self.branch_node] += veq;
        } else {
            // DC Operating Point: Inductor is a short circuit (0V drop)
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
        }
        Ok(())
    }
}
