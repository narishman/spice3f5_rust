use spice_core::circuit::Circuit;
use spice_core::devices::DeviceEnum;
use spice_core::devices::capacitor::Capacitor;
use spice_core::devices::diode::Diode;
use spice_core::devices::isource::CurrentSource;
use spice_core::devices::resistor::Resistor;
use spice_core::ni::do_tran;
use spice_core::node::NodeType;

#[test]
fn test_diode_transient_convergence() {
    // We recreate a non-linear low-pass circuit
    // A 10mA Current source charging a 1uF Capacitor in parallel with a Diode and 1k Resistor
    // Over time, voltage across node "out" builds until the diode exponentially clamps it.

    let mut circuit = Circuit::new();
    circuit.tran_max_iter = 100; // Increase transient iteration limits for exponential non-linear devices

    let node_out = circuit.get_or_create_node("out", NodeType::Voltage);
    let node_gnd = circuit.get_or_create_node("0", NodeType::Voltage);

    let mut devices = vec![
        DeviceEnum::CurrentSource(CurrentSource::new("I1", node_out, node_gnd, 0.01)), // 10 mA
        DeviceEnum::Resistor(Resistor::new("R1", node_out, node_gnd, 1000.0)),
        DeviceEnum::Capacitor(Capacitor::new("C1", node_out, node_gnd, 1e-6)), // 1 uF
        DeviceEnum::Diode(Diode::new("D1", node_out, node_gnd)), // Clamps roughly around 0.6V - 0.8V
    ];

    // Simulate for 1ms with 10us step sizes
    let step_size = 1e-5;
    let t_stop = 1e-3;

    let result = do_tran(&mut circuit, &mut devices, step_size, t_stop);
    assert!(
        result.is_ok(),
        "Transient analysis halted with error: {:?}",
        result.err()
    );

    // Evaluate final clamped voltage
    let v_out = circuit.rhs[node_out];

    // The 10mA source across 1k resistor would normally hit 10V.
    // The Diode should securely clamp the steady state voltage significantly lower (e.g. 0.3 - 0.7 V)
    assert!(
        v_out < 1.0,
        "Diode failed to clamp the output voltage! Evaluated: {}",
        v_out
    );
    assert!(
        v_out > 0.1,
        "Diode clamped inappropriately or solver collapsed: {}",
        v_out
    );
}
