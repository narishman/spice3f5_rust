use spice_core::circuit::Circuit;
use spice_core::loader::load_from_ast;
use spice_core::ni::do_dc_op;
use spice_parse::parse_line;

#[test]
fn test_spice_parser_integration() {
    let spice_netlist = [
        "* Norton Equivalent Test via Parser",
        "I1 in 0 1m",
        "R1 in out 10k",
        "R2 out 0 10k",
    ];

    let mut parsed_stmts = Vec::new();
    for line in spice_netlist.iter() {
        let (rem, stmt) = parse_line(line).expect("Failed to parse SPICE line");
        assert_eq!(rem, "", "Leftover junk in the parser");
        parsed_stmts.push(stmt);
    }

    let mut circuit = Circuit::new();

    // Map Ground Node explicitly
    circuit.get_or_create_node("0", spice_core::node::NodeType::Voltage);

    let mut devices =
        load_from_ast(&mut circuit, parsed_stmts).expect("Failed to generate AST Device tree");

    let result = do_dc_op(&mut circuit, &mut devices);
    assert!(result.is_ok(), "DC converged improperly");

    // "in" should be node 1, "out" should be node 2 because of parse order.
    // 0 is mapped immediately at start to idx 0.
    // However, nodes are stored as a Vec. Let's retrieve idx via iterators.
    let in_idx = circuit
        .nodes
        .iter()
        .find(|n| n.name == "in")
        .unwrap()
        .number;
    let out_idx = circuit
        .nodes
        .iter()
        .find(|n| n.name == "out")
        .unwrap()
        .number;

    let v_in = circuit.rhs[in_idx];
    let v_out = circuit.rhs[out_idx];

    assert!((v_in - 20.0).abs() < 1e-6);
    assert!((v_out - 10.0).abs() < 1e-6);
}
