use spice_core::circuit::Circuit;
use spice_core::devices::{
    DeviceEnum,
    bjt::{Bjt, BjtType},
    resistor::Resistor,
    vsource::VoltageSource,
};
use spice_core::ni::do_dc_op;

#[test]
fn test_bjt_dc_operating_point() {
    let mut circuit = Circuit::new();

    let base_node = circuit.get_or_create_node("B", spice_core::node::NodeType::Voltage);
    let col_node = circuit.get_or_create_node("C", spice_core::node::NodeType::Voltage);
    let emit_node = circuit.get_or_create_node("E", spice_core::node::NodeType::Voltage);
    let in_node = circuit.get_or_create_node("IN", spice_core::node::NodeType::Voltage);
    let vcc_node = circuit.get_or_create_node("VCC", spice_core::node::NodeType::Voltage);

    // In SPICE, ground is always node 0
    let gnd = 0;

    let mut devices = vec![
        // Voltage sources
        DeviceEnum::VoltageSource(VoltageSource::new("VCC", vcc_node, gnd, 5.0)), // 5V supply
        DeviceEnum::VoltageSource(VoltageSource::new("VIN", in_node, gnd, 1.5)),  // 1.5V input
        // Base Resistor: 100k
        DeviceEnum::Resistor(Resistor::new("RB", in_node, base_node, 100_000.0)),
        // Collector Resistor: 1k
        DeviceEnum::Resistor(Resistor::new("RC", vcc_node, col_node, 1_000.0)),
        // Emitter to ground directly
        DeviceEnum::VoltageSource(VoltageSource::new("VEE", emit_node, gnd, 0.0)),
    ];

    let mut bjt = Bjt::new("Q1", BjtType::NPN, col_node, base_node, emit_node, gnd);
    bjt.is_sat = 1e-16;
    bjt.bf = 100.0;

    devices.push(DeviceEnum::Bjt(bjt));

    do_dc_op(&mut circuit, &mut devices).expect("DC Operating Point failed to converge");

    let v_b = circuit.rhs[base_node];
    let v_c = circuit.rhs[col_node];
    let v_e = circuit.rhs[emit_node];

    println!("V_B: {}", v_b);
    println!("V_C: {}", v_c);
    println!("V_E: {}", v_e);

    // Check Forward Active Mode
    // V_BE ≈ 0.6 to 0.7 V
    assert!(
        v_b - v_e > 0.6 && v_b - v_e < 0.8,
        "V_BE should be around 0.7V, got {}",
        v_b - v_e
    );

    // V_BC < 0 (Reverse biased base-collector junction)
    assert!(
        v_b - v_c < 0.0,
        "Transistor should be in forward-active, but V_BC is {}",
        v_b - v_c
    );

    // Check beta amplification
    let i_b = (1.5 - v_b) / 100_000.0;
    let i_c = (5.0 - v_c) / 1_000.0;
    let beta = i_c / i_b;
    println!("Beta: {}", beta);
    assert!(
        (beta - 100.0).abs() < 5.0,
        "Beta should be close to 100, got {}",
        beta
    );
}
