use spice_core::circuit::Circuit;
use spice_core::devices::DeviceEnum;
use spice_core::devices::isource::CurrentSource;
use spice_core::devices::resistor::Resistor;
use spice_core::ni::do_dc_op;
use spice_core::node::NodeType;

#[test]
fn test_norton_equivalent_dc_op() {
    // We recreate a simple current-driven divider using Resistors and an I-Source:
    // Current Source (I1) = 0.001A (1mA) pushing into Node "in" from Ground "0"
    // R1 = 10k Ohm between "in" and "out"
    // R2 = 10k Ohm between "out" and "0"
    //
    // By Nodal Analysis (since R1 and R2 are in series, total R = 20k):
    // V_in = 1mA * 20k = 20.0V
    // V_out = 1mA * 10k = 10.0V

    let mut circuit = Circuit::new();

    let node_in = circuit.get_or_create_node("in", NodeType::Voltage);
    let node_out = circuit.get_or_create_node("out", NodeType::Voltage);
    let node_gnd = circuit.get_or_create_node("0", NodeType::Voltage);

    assert_eq!(
        node_gnd, 0,
        "Ground node must be mapped to internal matrix index 0"
    );

    let mut devices = vec![
        DeviceEnum::CurrentSource(CurrentSource::new("I1", node_in, node_gnd, 0.001)),
        DeviceEnum::Resistor(Resistor::new("R1", node_in, node_out, 10000.0)),
        DeviceEnum::Resistor(Resistor::new("R2", node_out, node_gnd, 10000.0)),
    ];

    let result = do_dc_op(&mut circuit, &mut devices);
    assert!(
        result.is_ok(),
        "DC Operating Point analysis halted with error: {:?}",
        result.err()
    );

    // Evaluate voltages
    let v_in = circuit.rhs[node_in];
    let v_out = circuit.rhs[node_out];

    // Assert expected limits over float ranges
    assert!(
        (v_in - 20.0).abs() < 1e-6,
        "V_in should be exactly 20.0V, got {}",
        v_in
    );
    assert!(
        (v_out - 10.0).abs() < 1e-6,
        "V_out should be exactly 10.0V, got {}",
        v_out
    );
}
