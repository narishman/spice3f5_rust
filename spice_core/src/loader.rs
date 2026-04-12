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
use spice_parse::ast::{DeviceStmt, Stmt};

/// Build a list of DeviceEnum models and populate Circuit Nodes from parsed ASTs
pub fn load_from_ast(circuit: &mut Circuit, stmts: Vec<Stmt>) -> Result<Vec<DeviceEnum>, String> {
    let mut devices = Vec::new();

    for stmt in stmts {
        match stmt {
            Stmt::Device(dev) => devices.push(load_device(circuit, dev)?),
            Stmt::Comment(_) => {}
            _ => {
                return Err(
                    "Control and Model statements are not yet supported for loading".to_string(),
                );
            }
        }
    }

    Ok(devices)
}

fn load_device(circuit: &mut Circuit, stmt: DeviceStmt) -> Result<DeviceEnum, String> {
    let ident = stmt.ident.to_uppercase();
    if ident.starts_with("R") {
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
    } else if ident.starts_with("C") {
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
    } else if ident.starts_with("L") {
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
    } else if ident.starts_with("V") {
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
    } else if ident.starts_with("I") {
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
    } else if ident.starts_with("D") {
        if stmt.nodes.len() != 2 {
            return Err(format!("Diode {} needs 2 nodes", stmt.ident));
        }
        let n1 = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let n2 = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        Ok(DeviceEnum::Diode(Diode::new(&stmt.ident, n1, n2)))
    } else if ident.starts_with("Q") {
        if stmt.nodes.len() != 3 {
            return Err(format!("BJT {} needs 3 nodes (C, B, E)", stmt.ident));
        }
        let c = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let b = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        let e = circuit.get_or_create_node(&stmt.nodes[2], NodeType::Voltage);
        Ok(DeviceEnum::Bjt(Bjt::new(
            &stmt.ident,
            BjtType::NPN,
            c,
            b,
            e,
            0,
        ))) // Ground assumed 0 for substrate right now
    } else if ident.starts_with("M") {
        if stmt.nodes.len() != 4 {
            return Err(format!("MOSFET {} needs 4 nodes (D, G, S, B)", stmt.ident));
        }
        let d = circuit.get_or_create_node(&stmt.nodes[0], NodeType::Voltage);
        let g = circuit.get_or_create_node(&stmt.nodes[1], NodeType::Voltage);
        let s = circuit.get_or_create_node(&stmt.nodes[2], NodeType::Voltage);
        let b = circuit.get_or_create_node(&stmt.nodes[3], NodeType::Voltage);
        Ok(DeviceEnum::Mos1(Mos1::new(
            &stmt.ident,
            MosType::NMOS,
            d,
            g,
            s,
            b,
        )))
    } else {
        Err(format!("Unsupported device ident: {}", stmt.ident))
    }
}
