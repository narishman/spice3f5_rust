use crate::circuit::Circuit;
use crate::devices::Device;
use crate::devices::devsup::dev_pnjlim;

const VT: f64 = 0.02585; // Thermal voltage at 300K

#[derive(PartialEq)]
pub enum BjtType {
    NPN,
    PNP,
}

pub struct Bjt {
    pub name: String,
    pub bjt_type: BjtType,

    // External Nodes
    pub col_node: usize,
    pub base_node: usize,
    pub emit_node: usize,
    pub subst_node: usize,

    // Internal Nodes
    pub col_prime_node: usize,
    pub base_prime_node: usize,
    pub emit_prime_node: usize,

    // Parameters
    pub area: f64,
    pub temp: f64,
    pub is_sat: f64,
    pub bf: f64,
    pub br: f64,
    pub ne: f64,
    pub nc: f64,
    pub vaf: f64,
    pub var: f64,
    pub ikf: f64,
    pub ikr: f64,
    pub ise: f64,
    pub isc: f64,
    pub rb: f64,
    pub rbm: f64,
    pub irb: f64,
    pub rc: f64,
    pub re: f64,

    // Derived vars
    pub state_idx: usize,

    // Matrix Pointers
    pub col_col_ptr: Option<usize>,
    pub base_base_ptr: Option<usize>,
    pub emit_emit_ptr: Option<usize>,
    pub col_prime_col_prime_ptr: Option<usize>,
    pub base_prime_base_prime_ptr: Option<usize>,
    pub emit_prime_emit_prime_ptr: Option<usize>,
    pub col_col_prime_ptr: Option<usize>,
    pub base_base_prime_ptr: Option<usize>,
    pub emit_emit_prime_ptr: Option<usize>,
    pub col_prime_col_ptr: Option<usize>,
    pub col_prime_base_prime_ptr: Option<usize>,
    pub col_prime_emit_prime_ptr: Option<usize>,
    pub base_prime_base_ptr: Option<usize>,
    pub base_prime_col_prime_ptr: Option<usize>,
    pub base_prime_emit_prime_ptr: Option<usize>,
    pub emit_prime_emit_ptr: Option<usize>,
    pub emit_prime_col_prime_ptr: Option<usize>,
    pub emit_prime_base_prime_ptr: Option<usize>,
    pub subst_subst_ptr: Option<usize>,
    pub col_prime_subst_ptr: Option<usize>,
    pub subst_col_prime_ptr: Option<usize>,
    pub base_col_prime_ptr: Option<usize>,
    pub col_prime_base_ptr: Option<usize>,
}

impl Bjt {
    pub fn new(
        name: impl Into<String>,
        bjt_type: BjtType,
        col_node: usize,
        base_node: usize,
        emit_node: usize,
        subst_node: usize,
    ) -> Self {
        Self {
            name: name.into(),
            bjt_type,
            col_node,
            base_node,
            emit_node,
            subst_node,

            // These will be allocated during setup
            col_prime_node: 0,
            base_prime_node: 0,
            emit_prime_node: 0,

            // Default model parameters
            area: 1.0,
            temp: 300.0,
            is_sat: 1e-16,
            bf: 100.0,
            br: 1.0,
            ne: 1.5,
            nc: 2.0,
            vaf: 0.0, // 0 means infinity (disabled)
            var: 0.0,
            ikf: 0.0,
            ikr: 0.0,
            ise: 0.0,
            isc: 0.0,
            rb: 0.0,
            rbm: 0.0,
            irb: 0.0,
            rc: 0.0,
            re: 0.0,

            state_idx: 0,
            col_col_ptr: None,
            base_base_ptr: None,
            emit_emit_ptr: None,
            col_prime_col_prime_ptr: None,
            base_prime_base_prime_ptr: None,
            emit_prime_emit_prime_ptr: None,
            col_col_prime_ptr: None,
            base_base_prime_ptr: None,
            emit_emit_prime_ptr: None,
            col_prime_col_ptr: None,
            col_prime_base_prime_ptr: None,
            col_prime_emit_prime_ptr: None,
            base_prime_base_ptr: None,
            base_prime_col_prime_ptr: None,
            base_prime_emit_prime_ptr: None,
            emit_prime_emit_ptr: None,
            emit_prime_col_prime_ptr: None,
            emit_prime_base_prime_ptr: None,
            subst_subst_ptr: None,
            col_prime_subst_ptr: None,
            subst_col_prime_ptr: None,
            base_col_prime_ptr: None,
            col_prime_base_ptr: None,
        }
    }
}

impl Device for Bjt {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        // Create internal prime nodes if resistances exist
        if self.rc > 0.0 {
            self.col_prime_node = circuit.get_or_create_node(
                &format!("{}_col_prime", self.name),
                crate::node::NodeType::Voltage,
            );
        } else {
            self.col_prime_node = self.col_node;
        }

        if self.rb > 0.0 {
            self.base_prime_node = circuit.get_or_create_node(
                &format!("{}_base_prime", self.name),
                crate::node::NodeType::Voltage,
            );
        } else {
            self.base_prime_node = self.base_node;
        }

        if self.re > 0.0 {
            self.emit_prime_node = circuit.get_or_create_node(
                &format!("{}_emit_prime", self.name),
                crate::node::NodeType::Voltage,
            );
        } else {
            self.emit_prime_node = self.emit_node;
        }

        self.state_idx = circuit.num_states;
        circuit.num_states += 21; // BJT uses 21 states internally per BJTload.c

        // Register all interacting sparse combinations
        self.col_col_ptr = circuit.matrix.get_element(self.col_node, self.col_node);
        self.base_base_ptr = circuit.matrix.get_element(self.base_node, self.base_node);
        self.emit_emit_ptr = circuit.matrix.get_element(self.emit_node, self.emit_node);
        self.col_prime_col_prime_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.col_prime_node);
        self.base_prime_base_prime_ptr = circuit
            .matrix
            .get_element(self.base_prime_node, self.base_prime_node);
        self.emit_prime_emit_prime_ptr = circuit
            .matrix
            .get_element(self.emit_prime_node, self.emit_prime_node);

        self.col_col_prime_ptr = circuit
            .matrix
            .get_element(self.col_node, self.col_prime_node);
        self.base_base_prime_ptr = circuit
            .matrix
            .get_element(self.base_node, self.base_prime_node);
        self.emit_emit_prime_ptr = circuit
            .matrix
            .get_element(self.emit_node, self.emit_prime_node);

        self.col_prime_col_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.col_node);
        self.col_prime_base_prime_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.base_prime_node);
        self.col_prime_emit_prime_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.emit_prime_node);

        self.base_prime_base_ptr = circuit
            .matrix
            .get_element(self.base_prime_node, self.base_node);
        self.base_prime_col_prime_ptr = circuit
            .matrix
            .get_element(self.base_prime_node, self.col_prime_node);
        self.base_prime_emit_prime_ptr = circuit
            .matrix
            .get_element(self.base_prime_node, self.emit_prime_node);

        self.emit_prime_emit_ptr = circuit
            .matrix
            .get_element(self.emit_prime_node, self.emit_node);
        self.emit_prime_col_prime_ptr = circuit
            .matrix
            .get_element(self.emit_prime_node, self.col_prime_node);
        self.emit_prime_base_prime_ptr = circuit
            .matrix
            .get_element(self.emit_prime_node, self.base_prime_node);

        self.subst_subst_ptr = circuit.matrix.get_element(self.subst_node, self.subst_node);
        self.col_prime_subst_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.subst_node);
        self.subst_col_prime_ptr = circuit
            .matrix
            .get_element(self.subst_node, self.col_prime_node);
        self.base_col_prime_ptr = circuit
            .matrix
            .get_element(self.base_node, self.col_prime_node);
        self.col_prime_base_ptr = circuit
            .matrix
            .get_element(self.col_prime_node, self.base_node);

        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        let type_multiplier = if self.bjt_type == BjtType::NPN {
            1.0
        } else {
            -1.0
        };

        let v_base_prime = if self.base_prime_node != 0 {
            circuit.rhs_old[self.base_prime_node]
        } else {
            0.0
        };
        let v_emit_prime = if self.emit_prime_node != 0 {
            circuit.rhs_old[self.emit_prime_node]
        } else {
            0.0
        };
        let v_col_prime = if self.col_prime_node != 0 {
            circuit.rhs_old[self.col_prime_node]
        } else {
            0.0
        };

        // Linear guesses for NR
        let v_be_matrix = type_multiplier * (v_base_prime - v_emit_prime);
        let v_bc_matrix = type_multiplier * (v_base_prime - v_col_prime);

        // PNJLIM bounded states
        let v_be_old = circuit.states[0][self.state_idx]; // BJTvbe
        let v_bc_old = circuit.states[0][self.state_idx + 1]; // BJTvbc

        let vcrit = VT * (VT / (1.414 * self.is_sat)).ln();

        let mut icheck_be = false;
        let mut icheck_bc = false;
        let vbe = dev_pnjlim(v_be_matrix, v_be_old, VT, vcrit, &mut icheck_be);
        let vbc = dev_pnjlim(v_bc_matrix, v_bc_old, VT, vcrit, &mut icheck_bc);

        if icheck_be || icheck_bc {
            circuit.noncon += 1;
        }

        // Core Gummel-Poon Equivalents
        let vtn = VT;

        // Determine base-emitter current and conductance
        let (cbe, gbe) = if vbe > -5.0 * VT {
            let evbe = (vbe / VT).exp();
            let cbe = self.is_sat * (evbe - 1.0) + circuit.gmin * vbe;
            let gbe = self.is_sat * evbe / VT + circuit.gmin;
            (cbe, gbe)
        } else {
            let gbe = -self.is_sat / vbe + circuit.gmin;
            (gbe * vbe, gbe)
        };

        // Determine base-collector current and conductance
        let (cbc, gbc) = if vbc > -5.0 * VT {
            let evbc = (vbc / VT).exp();
            let cbc = self.is_sat * (evbc - 1.0) + circuit.gmin * vbc;
            let gbc = self.is_sat * evbc / VT + circuit.gmin;
            (cbc, gbc)
        } else {
            let gbc = -self.is_sat / vbc + circuit.gmin;
            (gbc * vbc, gbc)
        };

        // Base Charge mapping Q_b (Early Effect)
        let q1 = 1.0
            / (1.0
                - (if self.vaf > 0.0 { vbc / self.vaf } else { 0.0 })
                - (if self.var > 0.0 { vbe / self.var } else { 0.0 }));

        let qb = q1; // Omitting High-Level Injection (ikf, ikr) for immediate scope, matching SPICE basic qb

        let dqbdve = if self.var > 0.0 {
            q1 * qb / self.var
        } else {
            0.0
        };
        let dqbdvc = if self.vaf > 0.0 {
            q1 * qb / self.vaf
        } else {
            0.0
        };

        let cc = (cbe - cbc) / qb - cbc / self.br;
        let cb = cbe / self.bf + cbc / self.br;

        let gx = if self.rb > 0.0 {
            1.0 / (self.rbm + (self.rb - self.rbm) / qb)
        } else {
            0.0
        };
        let gpi = gbe / self.bf;
        let gmu = gbc / self.br;

        let go = (gbc + (cbe - cbc) * dqbdvc / qb) / qb;
        let gm = (gbe - (cbe - cbc) * dqbdve / qb) / qb - go;

        let gcpr = if self.rc > 0.0 { 1.0 / self.rc } else { 0.0 };
        let gepr = if self.re > 0.0 { 1.0 / self.re } else { 0.0 };

        // Output State saving
        circuit.states[0][self.state_idx] = vbe;
        circuit.states[0][self.state_idx + 1] = vbc;

        // Current Excitations Setup
        let ceqbe = type_multiplier * (cc + cb - vbe * (gm + go + gpi) + vbc * go);
        let ceqbc = type_multiplier * (-cc + vbe * (gm + go) - vbc * (gmu + go));

        // Load Excitation Vector (RHS)
        if self.col_prime_node != 0 {
            circuit.rhs[self.col_prime_node] += ceqbc;
        }
        if self.base_prime_node != 0 {
            circuit.rhs[self.base_prime_node] += -ceqbe - ceqbc;
        }
        if self.emit_prime_node != 0 {
            circuit.rhs[self.emit_prime_node] += ceqbe;
        }

        // Load Jacobian Matrix
        if let Some(idx) = self.col_col_ptr {
            circuit.matrix.elements_arena[idx].real += gcpr;
        }
        if let Some(idx) = self.base_base_ptr {
            circuit.matrix.elements_arena[idx].real += gx;
        }
        if let Some(idx) = self.emit_emit_ptr {
            circuit.matrix.elements_arena[idx].real += gepr;
        }

        if let Some(idx) = self.col_prime_col_prime_ptr {
            circuit.matrix.elements_arena[idx].real += gmu + go + gcpr;
        }
        if let Some(idx) = self.base_prime_base_prime_ptr {
            circuit.matrix.elements_arena[idx].real += gx + gpi + gmu;
        }
        if let Some(idx) = self.emit_prime_emit_prime_ptr {
            circuit.matrix.elements_arena[idx].real += gpi + gepr + gm + go;
        }

        if let Some(idx) = self.col_col_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= gcpr;
        }
        if let Some(idx) = self.base_base_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= gx;
        }
        if let Some(idx) = self.emit_emit_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= gepr;
        }

        if let Some(idx) = self.col_prime_col_ptr {
            circuit.matrix.elements_arena[idx].real -= gcpr;
        }
        if let Some(idx) = self.col_prime_base_prime_ptr {
            circuit.matrix.elements_arena[idx].real += -gmu + gm;
        }
        if let Some(idx) = self.col_prime_emit_prime_ptr {
            circuit.matrix.elements_arena[idx].real += -gm - go;
        }

        if let Some(idx) = self.base_prime_base_ptr {
            circuit.matrix.elements_arena[idx].real -= gx;
        }
        if let Some(idx) = self.base_prime_col_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= gmu;
        }
        if let Some(idx) = self.base_prime_emit_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= gpi;
        }

        if let Some(idx) = self.emit_prime_emit_ptr {
            circuit.matrix.elements_arena[idx].real -= gepr;
        }
        if let Some(idx) = self.emit_prime_col_prime_ptr {
            circuit.matrix.elements_arena[idx].real -= go;
        }
        if let Some(idx) = self.emit_prime_base_prime_ptr {
            circuit.matrix.elements_arena[idx].real += -gpi - gm;
        }

        Ok(())
    }
}
