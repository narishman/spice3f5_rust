use crate::circuit::Circuit;
use crate::devices::Device;
use crate::devices::devsup::{dev_fetlim, dev_limvds, dev_pnjlim, dev_qmeyer};
use crate::integrate::integrate;

const VT: f64 = 0.02585; // Thermal voltage at 300K

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum MosType {
    NMOS,
    PMOS,
}

pub struct Mos1 {
    pub name: String,
    pub mos_type: MosType,

    // External Nodes
    pub d_node: usize,
    pub g_node: usize,
    pub s_node: usize,
    pub b_node: usize,

    // Internal Nodes
    pub d_node_prime: usize,
    pub s_node_prime: usize,

    // Parameters
    pub temp: f64,
    pub w: f64,
    pub l: f64,
    pub as_area: f64,
    pub ad_area: f64,

    // Device Conductances
    pub source_conductance: f64,
    pub drain_conductance: f64,

    // Internal State Tracker
    pub state_idx: usize,
    pub mode: i32, // 1 = normal, -1 = inverse
    pub off: bool,

    // Model Parameters (Tied to instance for direct calculation simplified)
    pub type_multiplier: f64,
    pub lat_diff: f64,
    pub t_sat_cur: f64,
    pub t_sat_cur_dens: f64,
    pub t_phi: f64,
    pub t_vto: f64,
    pub t_transconductance: f64,
    pub t_bulk_pot: f64,
    pub t_dep_cap: f64,
    pub t_vbi: f64,

    pub lambda: f64,
    pub gamma: f64,

    pub gate_source_overlap_cap_factor: f64,
    pub gate_drain_overlap_cap_factor: f64,
    pub gate_bulk_overlap_cap_factor: f64,
    pub oxide_cap_factor: f64,

    pub v_on: f64,
    pub v_dsat: f64,
    pub source_vcrit: f64,
    pub drain_vcrit: f64,

    pub cd: f64,
    pub cbs: f64,
    pub cbd: f64,
    pub gmbs: f64,
    pub gm: f64,
    pub gds: f64,
    pub gbd: f64,
    pub gbs: f64,

    pub capbd: f64,
    pub capbs: f64,
    pub c_bd: f64,
    pub c_bdsw: f64,
    pub c_bs: f64,
    pub c_bssw: f64,
    pub f2d: f64,
    pub f3d: f64,
    pub f4d: f64,
    pub f2s: f64,
    pub f3s: f64,
    pub f4s: f64,
    pub bulk_jct_bot_grading_coeff: f64,
    pub bulk_jct_side_grading_coeff: f64,

    pub capgs: f64,
    pub capgd: f64,
    pub capgb: f64,

    // Sparse Matrix Element Pointers
    pub dd_ptr: Option<usize>,
    pub gg_ptr: Option<usize>,
    pub ss_ptr: Option<usize>,
    pub bb_ptr: Option<usize>,
    pub dpdp_ptr: Option<usize>,
    pub spsp_ptr: Option<usize>,
    pub ddp_ptr: Option<usize>,
    pub gb_ptr: Option<usize>,
    pub gdp_ptr: Option<usize>,
    pub gsp_ptr: Option<usize>,
    pub ssp_ptr: Option<usize>,
    pub bdp_ptr: Option<usize>,
    pub bsp_ptr: Option<usize>,
    pub dpsp_ptr: Option<usize>,
    pub dpd_ptr: Option<usize>,
    pub bg_ptr: Option<usize>,
    pub dpg_ptr: Option<usize>,
    pub spg_ptr: Option<usize>,
    pub sps_ptr: Option<usize>,
    pub dpb_ptr: Option<usize>,
    pub spb_ptr: Option<usize>,
    pub spdp_ptr: Option<usize>,
}

impl Mos1 {
    pub fn new(
        name: impl Into<String>,
        mos_type: MosType,
        d_node: usize,
        g_node: usize,
        s_node: usize,
        b_node: usize,
    ) -> Self {
        Self {
            name: name.into(),
            mos_type,
            d_node,
            g_node,
            s_node,
            b_node,

            d_node_prime: 0,
            s_node_prime: 0,

            temp: 300.0,
            w: 1e-4,
            l: 1e-4,
            as_area: 0.0,
            ad_area: 0.0,

            source_conductance: 0.0,
            drain_conductance: 0.0,

            state_idx: 0,
            mode: 1,
            off: false,

            type_multiplier: if mos_type == MosType::NMOS { 1.0 } else { -1.0 },

            // Simplistic default parameters, mirroring basic MOS1 models
            lat_diff: 0.0,
            t_sat_cur: 1e-14,
            t_sat_cur_dens: 0.0,
            t_phi: 0.6,
            t_vto: if mos_type == MosType::NMOS { 1.0 } else { -1.0 },
            t_transconductance: 2.0e-5,
            t_bulk_pot: 0.8,
            t_dep_cap: 0.0,
            t_vbi: if mos_type == MosType::NMOS { 1.0 } else { -1.0 },

            lambda: 0.0,
            gamma: 0.0,
            gate_source_overlap_cap_factor: 0.0,
            gate_drain_overlap_cap_factor: 0.0,
            gate_bulk_overlap_cap_factor: 0.0,
            oxide_cap_factor: 0.0,

            v_on: 0.0,
            v_dsat: 0.0,
            source_vcrit: 0.0,
            drain_vcrit: 0.0,

            cd: 0.0,
            cbs: 0.0,
            cbd: 0.0,
            gmbs: 0.0,
            gm: 0.0,
            gds: 0.0,
            gbd: 0.0,
            gbs: 0.0,

            capbd: 0.0,
            capbs: 0.0,
            c_bd: 0.0,
            c_bdsw: 0.0,
            c_bs: 0.0,
            c_bssw: 0.0,
            f2d: 0.0,
            f3d: 0.0,
            f4d: 0.0,
            f2s: 0.0,
            f3s: 0.0,
            f4s: 0.0,
            bulk_jct_bot_grading_coeff: 0.5,
            bulk_jct_side_grading_coeff: 0.33,

            capgs: 0.0,
            capgd: 0.0,
            capgb: 0.0,

            dd_ptr: None,
            gg_ptr: None,
            ss_ptr: None,
            bb_ptr: None,
            dpdp_ptr: None,
            spsp_ptr: None,
            ddp_ptr: None,
            gb_ptr: None,
            gdp_ptr: None,
            gsp_ptr: None,
            ssp_ptr: None,
            bdp_ptr: None,
            bsp_ptr: None,
            dpsp_ptr: None,
            dpd_ptr: None,
            bg_ptr: None,
            dpg_ptr: None,
            spg_ptr: None,
            sps_ptr: None,
            dpb_ptr: None,
            spb_ptr: None,
            spdp_ptr: None,
        }
    }
}

impl Device for Mos1 {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        if self.drain_conductance > 0.0 {
            self.d_node_prime = circuit.get_or_create_node(
                &format!("{}_drain_prime", self.name),
                crate::node::NodeType::Voltage,
            );
        } else {
            self.d_node_prime = self.d_node;
        }

        if self.source_conductance > 0.0 {
            self.s_node_prime = circuit.get_or_create_node(
                &format!("{}_source_prime", self.name),
                crate::node::NodeType::Voltage,
            );
        } else {
            self.s_node_prime = self.s_node;
        }

        self.state_idx = circuit.num_states;
        circuit.num_states += 17;

        // Fetch matrix pointers
        self.dd_ptr = circuit.matrix.get_element(self.d_node, self.d_node);
        self.gg_ptr = circuit.matrix.get_element(self.g_node, self.g_node);
        self.ss_ptr = circuit.matrix.get_element(self.s_node, self.s_node);
        self.bb_ptr = circuit.matrix.get_element(self.b_node, self.b_node);

        self.dpdp_ptr = circuit
            .matrix
            .get_element(self.d_node_prime, self.d_node_prime);
        self.spsp_ptr = circuit
            .matrix
            .get_element(self.s_node_prime, self.s_node_prime);

        self.ddp_ptr = circuit.matrix.get_element(self.d_node, self.d_node_prime);
        self.gb_ptr = circuit.matrix.get_element(self.g_node, self.b_node);
        self.gdp_ptr = circuit.matrix.get_element(self.g_node, self.d_node_prime);
        self.gsp_ptr = circuit.matrix.get_element(self.g_node, self.s_node_prime);

        self.ssp_ptr = circuit.matrix.get_element(self.s_node, self.s_node_prime);
        self.bdp_ptr = circuit.matrix.get_element(self.b_node, self.d_node_prime);
        self.bsp_ptr = circuit.matrix.get_element(self.b_node, self.s_node_prime);

        self.dpsp_ptr = circuit
            .matrix
            .get_element(self.d_node_prime, self.s_node_prime);
        self.dpd_ptr = circuit.matrix.get_element(self.d_node_prime, self.d_node);
        self.bg_ptr = circuit.matrix.get_element(self.b_node, self.g_node);

        self.dpg_ptr = circuit.matrix.get_element(self.d_node_prime, self.g_node);
        self.spg_ptr = circuit.matrix.get_element(self.s_node_prime, self.g_node);

        self.sps_ptr = circuit.matrix.get_element(self.s_node_prime, self.s_node);
        self.dpb_ptr = circuit.matrix.get_element(self.d_node_prime, self.b_node);
        self.spb_ptr = circuit.matrix.get_element(self.s_node_prime, self.b_node);
        self.spdp_ptr = circuit
            .matrix
            .get_element(self.s_node_prime, self.d_node_prime);

        let vcrit = VT * (VT / (1.414 * self.t_sat_cur)).ln();
        self.source_vcrit = vcrit;
        self.drain_vcrit = vcrit;

        // Precompute flat-band / built-in threshold voltage at nominal temperature (1:1 with MOS1temp)
        self.type_multiplier = if self.mos_type == MosType::NMOS {
            1.0
        } else {
            -1.0
        };
        self.t_vbi = self.t_vto - self.type_multiplier * self.gamma * self.t_phi.sqrt();

        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        let vt = VT;
        let mut check = false;

        let effective_length = self.l - 2.0 * self.lat_diff;

        // SatCur evaluation
        let drain_sat_cur =
            if self.t_sat_cur_dens == 0.0 || self.ad_area == 0.0 || self.as_area == 0.0 {
                self.t_sat_cur
            } else {
                self.t_sat_cur_dens * self.ad_area
            };

        let source_sat_cur =
            if self.t_sat_cur_dens == 0.0 || self.ad_area == 0.0 || self.as_area == 0.0 {
                self.t_sat_cur
            } else {
                self.t_sat_cur_dens * self.as_area
            };

        let beta = self.t_transconductance * self.w / effective_length;

        let c_vbs = if self.b_node != 0 || self.s_node_prime != 0 {
            self.type_multiplier
                * (circuit.rhs_old[self.b_node] - circuit.rhs_old[self.s_node_prime])
        } else {
            0.0
        };

        let c_vgs = if self.g_node != 0 || self.s_node_prime != 0 {
            self.type_multiplier
                * (circuit.rhs_old[self.g_node] - circuit.rhs_old[self.s_node_prime])
        } else {
            0.0
        };

        let c_vds = if self.d_node_prime != 0 || self.s_node_prime != 0 {
            self.type_multiplier
                * (circuit.rhs_old[self.d_node_prime] - circuit.rhs_old[self.s_node_prime])
        } else {
            0.0
        };

        let mut vbs = c_vbs;
        let mut vgs = c_vgs;
        let mut vds = c_vds;
        let mut vbd = vbs - vds;
        let mut vgd = vgs - vds;

        let vbd_old = circuit.states[0][self.state_idx + 0]; // MOS1vbd
        let vbs_old = circuit.states[0][self.state_idx + 1]; // MOS1vbs
        let vgs_old = circuit.states[0][self.state_idx + 2]; // MOS1vgs
        let vds_old = circuit.states[0][self.state_idx + 3]; // MOS1vds
        let vgdo = vgs_old - vds_old;

        let von_prev = self.type_multiplier * self.v_on;

        // Node limiting to maintain numerical stability
        if vds_old >= 0.0 {
            vgs = dev_fetlim(vgs, vgs_old, von_prev);
            vds = vgs - vgd;
            vds = dev_limvds(vds, vds_old);
            vgd = vgs - vds;
        } else {
            vgd = dev_fetlim(vgd, vgdo, von_prev);
            vds = vgs - vgd;
            vds = -dev_limvds(-vds, -vds_old);
            vgs = vgd + vds;
        }

        if vds >= 0.0 {
            vbs = dev_pnjlim(vbs, vbs_old, vt, self.source_vcrit, &mut check);
            vbd = vbs - vds;
        } else {
            vbd = dev_pnjlim(vbd, vbd_old, vt, self.drain_vcrit, &mut check);
            vbs = vbd + vds;
        }

        if check {
            circuit.noncon += 1;
        }

        // Bulk-Source and Bulk-Drain diodes
        if vbs <= 0.0 {
            self.gbs = source_sat_cur / vt + circuit.gmin;
            self.cbs = self.gbs * vbs;
        } else {
            let evbs = (vbs / vt).min(40.0).exp();
            self.gbs = source_sat_cur * evbs / vt + circuit.gmin;
            self.cbs = source_sat_cur * (evbs - 1.0);
        }

        if vbd <= 0.0 {
            self.gbd = drain_sat_cur / vt + circuit.gmin;
            self.cbd = self.gbd * vbd;
        } else {
            let evbd = (vbd / vt).min(40.0).exp();
            self.gbd = drain_sat_cur * evbd / vt + circuit.gmin;
            self.cbd = drain_sat_cur * (evbd - 1.0);
        }

        self.mode = if vds >= 0.0 { 1 } else { -1 };

        let (cdrain, von, vdsat);
        let sarg;

        let bulk_v = if self.mode == 1 { vbs } else { vbd };

        if bulk_v <= 0.0 {
            sarg = (self.t_phi - bulk_v).max(0.0).sqrt();
        } else {
            let base_sarg = (self.t_phi).sqrt();
            sarg = (base_sarg - bulk_v / (base_sarg + base_sarg)).max(0.0);
        }

        von = (self.t_vbi * self.type_multiplier) + self.gamma * sarg;
        let vgst = (if self.mode == 1 { vgs } else { vgd }) - von;
        vdsat = vgst.max(0.0);

        let arg = if sarg <= 0.0 {
            0.0
        } else {
            self.gamma / (sarg + sarg)
        };

        if vgst <= 0.0 {
            // Cutoff
            cdrain = 0.0;
            self.gm = 0.0;
            self.gds = 0.0;
            self.gmbs = 0.0;
        } else {
            let betap = beta * (1.0 + self.lambda * (vds * (self.mode as f64)));
            if vgst <= (vds * (self.mode as f64)) {
                // Saturation
                cdrain = betap * vgst * vgst * 0.5;
                self.gm = betap * vgst;
                self.gds = self.lambda * beta * vgst * vgst * 0.5;
                self.gmbs = self.gm * arg;
            } else {
                // Linear
                cdrain =
                    betap * (vds * (self.mode as f64)) * (vgst - 0.5 * (vds * (self.mode as f64)));
                self.gm = betap * (vds * (self.mode as f64));
                self.gds = betap * (vgst - (vds * (self.mode as f64)))
                    + self.lambda
                        * beta
                        * (vds * (self.mode as f64))
                        * (vgst - 0.5 * (vds * (self.mode as f64)));
                self.gmbs = self.gm * arg;
            }
        }

        self.v_on = self.type_multiplier * von;
        self.v_dsat = self.type_multiplier * vdsat;
        self.cd = (self.mode as f64) * cdrain - self.cbd;

        // Store states for NR
        circuit.states[0][self.state_idx + 0] = vbd;
        circuit.states[0][self.state_idx + 1] = vbs;
        circuit.states[0][self.state_idx + 2] = vgs;
        circuit.states[0][self.state_idx + 3] = vds;

        // Meyer's capacitor model (1:1 with mos1load.c lines 821-955)
        let vgb = vgs - vbs;
        let oxide_cap = self.oxide_cap_factor * effective_length * self.w;
        let gs_ovl_cap = self.gate_source_overlap_cap_factor * self.w;
        let gd_ovl_cap = self.gate_drain_overlap_cap_factor * self.w;
        let gb_ovl_cap = self.gate_bulk_overlap_cap_factor * effective_length;

        let mut half_capgs = 0.0;
        let mut half_capgd = 0.0;
        let mut half_capgb = 0.0;
        if self.mode > 0 {
            dev_qmeyer(
                vgs,
                vgd,
                vgb,
                von,
                vdsat,
                &mut half_capgs,
                &mut half_capgd,
                &mut half_capgb,
                self.t_phi,
                oxide_cap,
            );
        } else {
            dev_qmeyer(
                vgd,
                vgs,
                vgb,
                von,
                vdsat,
                &mut half_capgd,
                &mut half_capgs,
                &mut half_capgb,
                self.t_phi,
                oxide_cap,
            );
        }
        circuit.states[0][self.state_idx + 4] = half_capgs;
        circuit.states[0][self.state_idx + 7] = half_capgd;
        circuit.states[0][self.state_idx + 10] = half_capgb;

        let mut gcgs = 0.0;
        let mut ceqgs = 0.0;
        let mut gcgd = 0.0;
        let mut ceqgd = 0.0;
        let mut gcgb = 0.0;
        let mut ceqgb = 0.0;

        if circuit.delta == 0.0 {
            let capgs = 2.0 * half_capgs + gs_ovl_cap;
            let capgd = 2.0 * half_capgd + gd_ovl_cap;
            let capgb = 2.0 * half_capgb + gb_ovl_cap;
            self.capgs = capgs;
            self.capgd = capgd;
            self.capgb = capgb;
            circuit.states[0][self.state_idx + 5] = vgs * capgs;
            circuit.states[0][self.state_idx + 8] = vgd * capgd;
            circuit.states[0][self.state_idx + 11] = vgb * capgb;
        } else {
            let vgs1 = circuit.states[1][self.state_idx + 2];
            let vgd1 = vgs1 - circuit.states[1][self.state_idx + 3];
            let vgb1 = vgs1 - circuit.states[1][self.state_idx + 1];

            let capgs = half_capgs + circuit.states[1][self.state_idx + 4] + gs_ovl_cap;
            let capgd = half_capgd + circuit.states[1][self.state_idx + 7] + gd_ovl_cap;
            let capgb = half_capgb + circuit.states[1][self.state_idx + 10] + gb_ovl_cap;
            self.capgs = capgs;
            self.capgd = capgd;
            self.capgb = capgb;

            circuit.states[0][self.state_idx + 5] =
                (vgs - vgs1) * capgs + circuit.states[1][self.state_idx + 5];
            circuit.states[0][self.state_idx + 8] =
                (vgd - vgd1) * capgd + circuit.states[1][self.state_idx + 8];
            circuit.states[0][self.state_idx + 11] =
                (vgb - vgb1) * capgb + circuit.states[1][self.state_idx + 11];

            if capgs > 0.0 {
                integrate(circuit, &mut gcgs, &mut ceqgs, capgs, self.state_idx + 5);
                let ag0 = 2.0 / circuit.delta;
                ceqgs = ceqgs - gcgs * vgs + ag0 * circuit.states[0][self.state_idx + 5];
            } else {
                circuit.states[0][self.state_idx + 6] = 0.0;
            }
            if capgd > 0.0 {
                integrate(circuit, &mut gcgd, &mut ceqgd, capgd, self.state_idx + 8);
                let ag0 = 2.0 / circuit.delta;
                ceqgd = ceqgd - gcgd * vgd + ag0 * circuit.states[0][self.state_idx + 8];
            } else {
                circuit.states[0][self.state_idx + 9] = 0.0;
            }
            if capgb > 0.0 {
                integrate(circuit, &mut gcgb, &mut ceqgb, capgb, self.state_idx + 11);
                let ag0 = 2.0 / circuit.delta;
                ceqgb = ceqgb - gcgb * vgb + ag0 * circuit.states[0][self.state_idx + 11];
            } else {
                circuit.states[0][self.state_idx + 12] = 0.0;
            }
        }

        // Compute equivalents
        let ceqbs = self.type_multiplier * (self.cbs - (self.gbs - circuit.gmin) * vbs);
        let ceqbd = self.type_multiplier * (self.cbd - (self.gbd - circuit.gmin) * vbd);

        let (xnrm, xrev, cdreq);
        if self.mode >= 0 {
            xnrm = 1.0;
            xrev = 0.0;
            cdreq =
                self.type_multiplier * (cdrain - self.gds * vds - self.gm * vgs - self.gmbs * vbs);
        } else {
            xnrm = 0.0;
            xrev = 1.0;
            cdreq = -self.type_multiplier
                * (cdrain - self.gds * (-vds) - self.gm * vgd - self.gmbs * vbd);
        }

        // Load Excitation current vector (RHS)
        if self.g_node != 0 {
            circuit.rhs[self.g_node] -= self.type_multiplier * (ceqgs + ceqgb + ceqgd);
        }
        if self.b_node != 0 {
            circuit.rhs[self.b_node] -= ceqbs + ceqbd - self.type_multiplier * ceqgb;
        }
        if self.d_node_prime != 0 {
            circuit.rhs[self.d_node_prime] += ceqbd - cdreq + self.type_multiplier * ceqgd;
        }
        if self.s_node_prime != 0 {
            circuit.rhs[self.s_node_prime] += cdreq + ceqbs + self.type_multiplier * ceqgs;
        }

        // Load Jacobian Matrix (Y-Matrix)
        if let Some(idx) = self.dd_ptr {
            circuit.matrix.elements_arena[idx].real += self.drain_conductance;
        }
        if let Some(idx) = self.gg_ptr {
            circuit.matrix.elements_arena[idx].real += gcgd + gcgs + gcgb;
        }
        if let Some(idx) = self.ss_ptr {
            circuit.matrix.elements_arena[idx].real += self.source_conductance;
        }
        if let Some(idx) = self.bb_ptr {
            circuit.matrix.elements_arena[idx].real += self.gbd + self.gbs + gcgb;
        }

        if let Some(idx) = self.dpdp_ptr {
            circuit.matrix.elements_arena[idx].real +=
                self.drain_conductance + self.gds + self.gbd + xrev * (self.gm + self.gmbs) + gcgd;
        }
        if let Some(idx) = self.spsp_ptr {
            circuit.matrix.elements_arena[idx].real +=
                self.source_conductance + self.gds + self.gbs + xnrm * (self.gm + self.gmbs) + gcgs;
        }

        if let Some(idx) = self.ddp_ptr {
            circuit.matrix.elements_arena[idx].real -= self.drain_conductance;
        }
        if let Some(idx) = self.gb_ptr {
            circuit.matrix.elements_arena[idx].real -= gcgb;
        }
        if let Some(idx) = self.gdp_ptr {
            circuit.matrix.elements_arena[idx].real -= gcgd;
        }
        if let Some(idx) = self.gsp_ptr {
            circuit.matrix.elements_arena[idx].real -= gcgs;
        }
        if let Some(idx) = self.ssp_ptr {
            circuit.matrix.elements_arena[idx].real -= self.source_conductance;
        }

        if let Some(idx) = self.bg_ptr {
            circuit.matrix.elements_arena[idx].real -= gcgb;
        }
        if let Some(idx) = self.bdp_ptr {
            circuit.matrix.elements_arena[idx].real -= self.gbd;
        }
        if let Some(idx) = self.bsp_ptr {
            circuit.matrix.elements_arena[idx].real -= self.gbs;
        }
        if let Some(idx) = self.dpd_ptr {
            circuit.matrix.elements_arena[idx].real -= self.drain_conductance;
        }

        if let Some(idx) = self.dpg_ptr {
            circuit.matrix.elements_arena[idx].real += (xnrm - xrev) * self.gm - gcgd;
        }
        if let Some(idx) = self.dpb_ptr {
            circuit.matrix.elements_arena[idx].real += -self.gbd + (xnrm - xrev) * self.gmbs;
        }
        if let Some(idx) = self.dpsp_ptr {
            circuit.matrix.elements_arena[idx].real += -self.gds - xnrm * (self.gm + self.gmbs);
        }

        if let Some(idx) = self.spg_ptr {
            circuit.matrix.elements_arena[idx].real += -(xnrm - xrev) * self.gm - gcgs;
        }
        if let Some(idx) = self.sps_ptr {
            circuit.matrix.elements_arena[idx].real -= self.source_conductance;
        }
        if let Some(idx) = self.spb_ptr {
            circuit.matrix.elements_arena[idx].real += -self.gbs - (xnrm - xrev) * self.gmbs;
        }
        if let Some(idx) = self.spdp_ptr {
            circuit.matrix.elements_arena[idx].real += -self.gds - xrev * (self.gm + self.gmbs);
        }

        Ok(())
    }
}
