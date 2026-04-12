use crate::circuit::Circuit;
use crate::devices::Device;
use crate::devices::devsup::dev_pnjlim;

const VT: f64 = 0.02585; // Thermal voltage at 300K

pub struct Diode {
    pub name: String,
    pub pos_node: usize,
    pub neg_node: usize,

    // Model parameters
    pub saturation_current: f64, // IS
    pub emission_coeff: f64,     // N

    pub state_idx: usize, // Pointer to active solver state (for voltage bounding)
    pub pos_pos_ptr: Option<usize>,
    pub neg_neg_ptr: Option<usize>,
    pub pos_neg_ptr: Option<usize>,
    pub neg_pos_ptr: Option<usize>,
}

impl Diode {
    pub fn new(name: impl Into<String>, pos_node: usize, neg_node: usize) -> Self {
        Self {
            name: name.into(),
            pos_node,
            neg_node,
            saturation_current: 1e-14,
            emission_coeff: 1.0,
            state_idx: 0,
            pos_pos_ptr: None,
            neg_neg_ptr: None,
            pos_neg_ptr: None,
            neg_pos_ptr: None,
        }
    }
}

impl Device for Diode {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        self.state_idx = circuit.num_states;
        circuit.num_states += 1; // 1 state for Diode Voltage (Newton-Raphson tracking)

        self.pos_pos_ptr = circuit.matrix.get_element(self.pos_node, self.pos_node);
        self.neg_neg_ptr = circuit.matrix.get_element(self.neg_node, self.neg_node);
        self.pos_neg_ptr = circuit.matrix.get_element(self.pos_node, self.neg_node);
        self.neg_pos_ptr = circuit.matrix.get_element(self.neg_node, self.pos_node);
        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        let v_pos = if self.pos_node != 0 {
            circuit.rhs_old[self.pos_node]
        } else {
            0.0
        };
        let v_neg = if self.neg_node != 0 {
            circuit.rhs_old[self.neg_node]
        } else {
            0.0
        };

        // This is the linear solver guess for current iteration
        let v_d_matrix = v_pos - v_neg;

        // This is the accepted/clamped voltage from the PREVIOUS iteration
        let v_d_old = circuit.states[0][self.state_idx];

        let vt_n = VT * self.emission_coeff;

        // Compute critical voltage dynamically based on thermal factors and saturation current
        let vcrit = vt_n * (vt_n / (1.414 * self.saturation_current)).ln();

        let mut icheck = false;
        let v_d_new = dev_pnjlim(v_d_matrix, v_d_old, vt_n, vcrit, &mut icheck);

        // Output convergence status based on icheck
        if icheck {
            circuit.noncon += 1;
        }

        // I_D = I_S * (e^(V_D / (N * V_T)) - 1)
        let exp_term = (v_d_new / vt_n).exp();
        let i_d = self.saturation_current * (exp_term - 1.0);

        // Conductance G = dI_D / dV_D = I_S / (N * V_T) * e^(V_D / (N * V_T))
        let geq = (self.saturation_current / vt_n) * exp_term;

        // Equivalent current source for N-R: I_eq = I_D - G * V_D
        let ceq = i_d - geq * v_d_new;

        // Stamp Conductance G into MNA Matrix
        if let Some(idx) = self.pos_pos_ptr {
            circuit.matrix.elements_arena[idx].real += geq;
        }
        if let Some(idx) = self.neg_neg_ptr {
            circuit.matrix.elements_arena[idx].real += geq;
        }
        if let Some(idx) = self.pos_neg_ptr {
            circuit.matrix.elements_arena[idx].real -= geq;
        }
        if let Some(idx) = self.neg_pos_ptr {
            circuit.matrix.elements_arena[idx].real -= geq;
        }

        // Stamp Equivalent Current I_eq into RHS Vector
        if self.pos_node != 0 {
            circuit.rhs[self.pos_node] -= ceq;
        }
        if self.neg_node != 0 {
            circuit.rhs[self.neg_node] += ceq;
        }

        // Save computed voltage for NEXT iteration mapping
        circuit.states[0][self.state_idx] = v_d_new;

        Ok(())
    }
}
