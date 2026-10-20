use spice_core::circuit::Circuit;
use spice_core::devices::{
    DeviceEnum,
    mos1::{Mos1, MosType},
    resistor::Resistor,
    vsource::VoltageSource,
};
use spice_core::ni::do_dc_op;
use spice_core::node::NodeType;

#[test]
fn test_nmos_saturation_dc_op() {
    // NMOS common-source stage:
    // VDD = 5V, RD = 10k between VDD and Drain
    // VG = 3V, Source = 0V, Bulk = 0V
    // Model defaults: VTO = 1.0V, KP = 2e-5 A/V^2, W = 100u, L = 100u -> Beta = 2e-5 A/V^2
    // Since VGS = 3V > VTO = 1V, Vgst = 2V.
    // Saturation current: Id = 0.5 * Beta * Vgst^2 = 0.5 * 2e-5 * 4 = 40 uA.
    // V_D = VDD - Id * RD = 5.0 - 40e-6 * 10e3 = 4.6V.
    // Since V_DS = 4.6V >= V_dsat (2.0V), NMOS is in saturation!
    let mut circuit = Circuit::new();

    let vdd = circuit.get_or_create_node("vdd", NodeType::Voltage);
    let drain = circuit.get_or_create_node("drain", NodeType::Voltage);
    let gate = circuit.get_or_create_node("gate", NodeType::Voltage);
    let gnd = 0;

    let mos = Mos1::new("M1", MosType::NMOS, drain, gate, gnd, gnd);

    let mut devices = vec![
        DeviceEnum::VoltageSource(VoltageSource::new("VDD", vdd, gnd, 5.0)),
        DeviceEnum::VoltageSource(VoltageSource::new("VG", gate, gnd, 3.0)),
        DeviceEnum::Resistor(Resistor::new("RD", vdd, drain, 10_000.0)),
        DeviceEnum::Mos1(mos),
    ];

    do_dc_op(&mut circuit, &mut devices).expect("NMOS DC OP failed to converge");

    let v_drain = circuit.rhs[drain];
    assert!(
        (v_drain - 4.6).abs() < 1e-3,
        "Expected V_drain ≈ 4.6V in saturation, got {}",
        v_drain
    );
}

#[test]
fn test_cmos_inverter_dc_op() {
    // CMOS Inverter with VDD = 5.0V:
    // Case 1: Vin = 0.0V -> Vout ≈ 5.0V
    // Case 2: Vin = 5.0V -> Vout ≈ 0.0V
    for &(v_in, expected_v_out) in &[(0.0, 5.0), (5.0, 0.0)] {
        let mut circuit = Circuit::new();
        let vdd = circuit.get_or_create_node("vdd", NodeType::Voltage);
        let in_node = circuit.get_or_create_node("in", NodeType::Voltage);
        let out_node = circuit.get_or_create_node("out", NodeType::Voltage);
        let gnd = 0;

        let nmos = Mos1::new("MN", MosType::NMOS, out_node, in_node, gnd, gnd);
        let pmos = Mos1::new("MP", MosType::PMOS, out_node, in_node, vdd, vdd);

        let mut devices = vec![
            DeviceEnum::VoltageSource(VoltageSource::new("VDD", vdd, gnd, 5.0)),
            DeviceEnum::VoltageSource(VoltageSource::new("VIN", in_node, gnd, v_in)),
            DeviceEnum::Mos1(nmos),
            DeviceEnum::Mos1(pmos),
        ];

        do_dc_op(&mut circuit, &mut devices).expect("CMOS Inverter DC OP failed to converge");

        let v_out = circuit.rhs[out_node];
        assert!(
            (v_out - expected_v_out).abs() < 1e-2,
            "For Vin={}, expected Vout ≈ {}, got {}",
            v_in,
            expected_v_out,
            v_out
        );
    }
}
