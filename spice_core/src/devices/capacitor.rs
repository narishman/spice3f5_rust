use crate::circuit::Circuit;
use crate::devices::Device;
use crate::integrate::integrate;

pub struct Capacitor {
    pub ident: String,
    pub n1: usize,
    pub n2: usize,
    pub capacitance: f64,
    pub state_idx: usize, // Pointer to state buffer slot

    // Matrix Elements
    pub m_n1_n1: Option<usize>,
    pub m_n2_n2: Option<usize>,
    pub m_n1_n2: Option<usize>,
    pub m_n2_n1: Option<usize>,
}

impl Capacitor {
    pub fn new(ident: &str, n1: usize, n2: usize, capacitance: f64) -> Self {
        Self {
            ident: ident.to_string(),
            n1,
            n2,
            capacitance,
            state_idx: 0,
            m_n1_n1: None,
            m_n2_n2: None,
            m_n1_n2: None,
            m_n2_n1: None,
        }
    }
}

impl Device for Capacitor {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        // Allocate 2 states for Capacitor: Charge (Q) and Current (I)
        self.state_idx = circuit.num_states;
        circuit.num_states += 2;

        if self.n1 != 0 {
            self.m_n1_n1 = circuit.matrix.get_element(self.n1, self.n1);
        }
        if self.n2 != 0 {
            self.m_n2_n2 = circuit.matrix.get_element(self.n2, self.n2);
        }
        if self.n1 != 0 && self.n2 != 0 {
            self.m_n1_n2 = circuit.matrix.get_element(self.n1, self.n2);
            self.m_n2_n1 = circuit.matrix.get_element(self.n2, self.n1);
        }

        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        // If transient step
        if circuit.delta > 0.0 {
            let v_n1 = if self.n1 != 0 {
                circuit.rhs_old[self.n1]
            } else {
                0.0
            };
            let v_n2 = if self.n2 != 0 {
                circuit.rhs_old[self.n2]
            } else {
                0.0
            };
            let v_cap = v_n1 - v_n2;

            // Current charge Q = C * V
            circuit.states[0][self.state_idx] = self.capacitance * v_cap;

            let mut geq = 0.0;
            let mut ceq = 0.0;

            integrate(
                circuit,
                &mut geq,
                &mut ceq,
                self.capacitance,
                self.state_idx,
            );

            // Stamp into MNA
            if let Some(idx) = self.m_n1_n1 {
                circuit.matrix.elements_arena[idx].real += geq;
            }
            if let Some(idx) = self.m_n2_n2 {
                circuit.matrix.elements_arena[idx].real += geq;
            }
            if let Some(idx) = self.m_n1_n2 {
                circuit.matrix.elements_arena[idx].real -= geq;
            }
            if let Some(idx) = self.m_n2_n1 {
                circuit.matrix.elements_arena[idx].real -= geq;
            }

            if self.n1 != 0 {
                circuit.rhs[self.n1] -= ceq;
            }
            if self.n2 != 0 {
                circuit.rhs[self.n2] += ceq;
            }
        }
        Ok(())
    }
}
