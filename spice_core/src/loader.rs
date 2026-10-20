use crate::circuit::Circuit;
use crate::devices::DeviceEnum;
use crate::devices::bjt::{Bjt, BjtType};
use crate::devices::capacitor::Capacitor;
use crate::devices::diode::Diode;
use crate::devices::inductor::Inductor;
use crate::devices::isource::CurrentSource;
use crate::devices::mos1::{Mos1, MosType};
use crate::devices::resistor::Resistor;
use crate::devices::vsource::VoltageSource;
use crate::node::NodeType;
use spice_parse::ast::{DeviceStmt, ModelStmt, Stmt};
use std::collections::HashMap;

/// Build a list of DeviceEnum models and populate Circuit Nodes from parsed ASTs
pub fn load_from_ast(circuit: &mut Circuit, stmts: Vec<Stmt>) -> Result<Vec<DeviceEnum>, String> {
    let mut models: HashMap<String, ModelStmt> = HashMap::new();
    let mut device_stmts: Vec<DeviceStmt> = Vec::new();

    for stmt in stmts {
        match stmt {
            Stmt::Device(dev) => device_stmts.push(dev),
            Stmt::Model(m) => {
                models.insert(m.name.to_uppercase(), m);
            }
            Stmt::Comment(_) | Stmt::Control(_) => {}
        }
    }

    let mut devices = Vec::new();
    for dev in device_stmts {
        devices.push(load_device(circuit, dev, &models)?);
    }

    Ok(devices)
}

fn load_device(
    circuit: &mut Circuit,
    stmt: DeviceStmt,
    models: &HashMap<String, ModelStmt>,
) -> Result<DeviceEnum, String> {
    let ident = stmt.ident.to_uppercase();
    let model_def = stmt
        .model
        .as_ref()
        .and_then(|m| models.get(&m.to_uppercase()));

    if ident.starts_with('R') {
        if stmt.nodes.len() != 2 {
            return Err(format!("Resistor {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::Resistor(Resistor::new(
            &stmt.ident,
            n1,
            n2,
            stmt.value,
        )))
    } else if ident.starts_with('C') {
        if stmt.nodes.len() != 2 {
            return Err(format!("Capacitor {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::Capacitor(Capacitor::new(
            &stmt.ident,
            n1,
            n2,
            stmt.value,
        )))
    } else if ident.starts_with('L') {
        if stmt.nodes.len() != 2 {
            return Err(format!("Inductor {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::Inductor(Inductor::new(
            &stmt.ident,
            n1,
            n2,
            stmt.value,
        )))
    } else if ident.starts_with('V') {
        if stmt.nodes.len() != 2 {
            return Err(format!("VSource {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::VoltageSource(VoltageSource::new(
            &stmt.ident,
            n1,
            n2,
            stmt.value,
        )))
    } else if ident.starts_with('I') {
        if stmt.nodes.len() != 2 {
            return Err(format!("ISource {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::CurrentSource(CurrentSource::new(
            &stmt.ident,
            n1,
            n2,
            stmt.value,
        )))
    } else if ident.starts_with('D') {
        if stmt.nodes.len() != 2 {
            return Err(format!("Diode {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        let mut diode = Diode::new(&stmt.ident, n1, n2);
        if let Some(m) = model_def {
            for (k, v) in &m.params {
                match k.as_str() {
                    "IS" => diode.saturation_current = *v,
                    "N" => diode.emission_coeff = *v,
                    _ => {}
                }
            }
        }
        for (k, v) in &stmt.params {
            match k.as_str() {
                "IS" => diode.saturation_current = *v,
                "N" => diode.emission_coeff = *v,
                _ => {}
            }
        }
        Ok(DeviceEnum::Diode(diode))
    } else if ident.starts_with('Q') {
        if stmt.nodes.len() != 3 {
            return Err(format!("BJT {} needs 3 nodes (C, B, E)", stmt.ident));
        }
        let c = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let b = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        let e = circuit.get_or_create_node(&stmt.nodes[2], NodeType::Voltage);
        let bjt_type = match model_def.map(|m| m.model_type.as_str()) {
            Some("PNP") => BjtType::PNP,
            _ => match stmt.model.as_deref().map(|s| s.to_uppercase()) {
                Some(ref s) if s == "PNP" => BjtType::PNP,
                _ => BjtType::NPN,
            },
        };
        let mut bjt = Bjt::new(&stmt.ident, bjt_type, c, b, e, 0);
        if let Some(m) = model_def {
            apply_bjt_params(&mut bjt, &m.params);
        }
        apply_bjt_params(&mut bjt, &stmt.params);
        Ok(DeviceEnum::Bjt(bjt))
    } else if ident.starts_with('M') {
        if stmt.nodes.len() != 4 {
            return Err(format!("MOSFET {} needs 4 nodes (D, G, S, B)", stmt.ident));
        }
        let d = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let g = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        let s = circuit.get_or_create_node(&stmt.nodes[2], NodeType::Voltage);
        let b = circuit.get_or_create_node(&stmt.nodes[3], NodeType::Voltage);
        let mos_type = match model_def.map(|m| m.model_type.as_str()) {
            Some("PMOS") => MosType::PMOS,
            _ => match stmt.model.as_deref().map(|s| s.to_uppercase()) {
                Some(ref s) if s == "PMOS" => MosType::PMOS,
                _ => MosType::NMOS,
            },
        };
        let mut mos = Mos1::new(&stmt.ident, mos_type, d, g, s, b);
        if let Some(m) = model_def {
            apply_mos1_params(&mut mos, &m.params);
        }
        apply_mos1_params(&mut mos, &stmt.params);
        Ok(DeviceEnum::Mos1(mos))
    } else {
        Err(format!("Unsupported device ident: {}", stmt.ident))
    }
}

fn apply_bjt_params(bjt: &mut Bjt, params: &[(String, f64)]) {
    for (k, v) in params {
        match k.as_str() {
            "IS" => bjt.is_sat = *v,
            "BF" => bjt.bf = *v,
            "BR" => bjt.br = *v,
            "NE" => bjt.ne = *v,
            "NC" => bjt.nc = *v,
            "VAF" | "VA" => bjt.vaf = *v,
            "VAR" => bjt.var = *v,
            "IKF" => bjt.ikf = *v,
            "IKR" => bjt.ikr = *v,
            "ISE" => bjt.ise = *v,
            "ISC" => bjt.isc = *v,
            "RB" => bjt.rb = *v,
            "RBM" => bjt.rbm = *v,
            "IRB" => bjt.irb = *v,
            "RC" => bjt.rc = *v,
            "RE" => bjt.re = *v,
            "AREA" => bjt.area = *v,
            _ => {}
        }
    }
}

fn apply_mos1_params(mos: &mut Mos1, params: &[(String, f64)]) {
    for (k, v) in params {
        match k.as_str() {
            "W" => mos.w = *v,
            "L" => mos.l = *v,
            "AS" => mos.as_area = *v,
            "AD" => mos.ad_area = *v,
            "VTO" | "VT0" => mos.t_vto = *v,
            "KP" => mos.t_transconductance = *v,
            "GAMMA" => mos.gamma = *v,
            "PHI" => mos.t_phi = *v,
            "LAMBDA" => mos.lambda = *v,
            "RD" => mos.drain_conductance = if *v > 0.0 { 1.0 / *v } else { 0.0 },
            "RS" => mos.source_conductance = if *v > 0.0 { 1.0 / *v } else { 0.0 },
            "IS" => mos.t_sat_cur = *v,
            "JS" => mos.t_sat_cur_dens = *v,
            "PB" => mos.t_bulk_pot = *v,
            "LD" => mos.lat_diff = *v,
            "CGSO" => mos.gate_source_overlap_cap_factor = *v,
            "CGDO" => mos.gate_drain_overlap_cap_factor = *v,
            "CGBO" => mos.gate_bulk_overlap_cap_factor = *v,
            "TOX" => {
                if *v > 0.0 {
                    mos.oxide_cap_factor = 3.9 * 8.854214871e-12 / *v;
                }
            }
            _ => {}
        }
    }
}
