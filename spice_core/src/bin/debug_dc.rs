use spice_core::circuit::Circuit;
use spice_core::devices::{
    Device, DeviceEnum,
    bjt::{Bjt, BjtType},
    resistor::Resistor,
    vsource::VoltageSource,
};
use spice_core::ni::do_dc_op;

fn main() {
    let mut circuit = Circuit::new();

    let base_node = circuit.get_or_create_node("B", spice_core::node::NodeType::Voltage);
    let col_node = circuit.get_or_create_node("C", spice_core::node::NodeType::Voltage);
    let emit_node = circuit.get_or_create_node("E", spice_core::node::NodeType::Voltage);
    let in_node = circuit.get_or_create_node("IN", spice_core::node::NodeType::Voltage);
    let vcc_node = circuit.get_or_create_node("VCC", spice_core::node::NodeType::Voltage);

    let gnd = 0;

    let mut devices = vec![
        DeviceEnum::VoltageSource(VoltageSource::new("VCC", vcc_node, gnd, 5.0)),
        DeviceEnum::VoltageSource(VoltageSource::new("VIN", in_node, gnd, 1.5)),
        DeviceEnum::Resistor(Resistor::new("RB", in_node, base_node, 100_000.0)),
        DeviceEnum::Resistor(Resistor::new("RC", vcc_node, col_node, 1_000.0)),
        DeviceEnum::VoltageSource(VoltageSource::new("VEE", emit_node, gnd, 0.0)),
    ];

    let mut bjt = Bjt::new("Q1", BjtType::NPN, col_node, base_node, emit_node, gnd);
    bjt.is_sat = 1e-16;
    bjt.bf = 100.0;
    devices.push(DeviceEnum::Bjt(bjt));

    for dev in devices.iter_mut() {
        dev.setup(&mut circuit).unwrap();
    }

    println!("Matrix Mapping after Setup:");
    for i in 1..=circuit.matrix.size {
        println!(
            "Internal Row {}: Ext Row {}, Internal Col {}: Ext Col {}",
            i, circuit.matrix.int_to_ext_row_map[i], i, circuit.matrix.int_to_ext_col_map[i]
        );
    }

    println!("Initial Matrix Elements BEFORE do_dc_op:");
    for i in 1..=circuit.matrix.size {
        let mut p = circuit.matrix.first_in_col[i];
        while let Some(idx) = p {
            let elem = &circuit.matrix.elements_arena[idx];
            println!(
                "IntRow {} IntCol {} -> real: {}",
                elem.row, elem.col, elem.real
            );
            p = elem.next_in_col;
        }
    }

    // Let ni do the work properly
    do_dc_op(&mut circuit, &mut devices).unwrap();

    println!("States output from ni::do_dc_op (circuit.rhs):");
    for (i, val) in circuit.rhs.iter().enumerate() {
        println!("{}: {}", i, val);
    }

    println!("Matrix Mapping after do_dc_op:");
    for i in 1..=circuit.matrix.size {
        println!(
            "Internal Row {}: Ext Row {}, Internal Col {}: Ext Col {}",
            i, circuit.matrix.int_to_ext_row_map[i], i, circuit.matrix.int_to_ext_col_map[i]
        );
    }

    println!("Factored Matrix Elements:");
    for i in 1..=circuit.matrix.size {
        let mut p = circuit.matrix.first_in_row[i];
        while let Some(idx) = p {
            let elem = &circuit.matrix.elements_arena[idx];
            println!(
                "IntRow {} IntCol {} -> Ext({} {}) -> real: {:?}",
                elem.row,
                elem.col,
                circuit.matrix.int_to_ext_row_map[elem.row],
                circuit.matrix.int_to_ext_col_map[elem.col],
                elem.real
            );
            p = elem.next_in_row;
        }
    }
}
