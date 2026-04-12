use spice_math::spdefs::MatrixFrame;
use spice_math::spfactor::PartitionMode;

const EPSILON: f64 = 1e-9;

/// Tests the basic allocation, element setup, and clearing operations of the sparse matrix engine.
#[test]
fn test_matrix_allocation_and_setup() {
    let mut matrix = MatrixFrame::new(0, false);

    // Test size scaling via translation maps (1:1 with Berkeley API paradigm)
    let idx1 = matrix.get_element(1, 1).unwrap();
    let idx2 = matrix.get_element(2, 2).unwrap();
    let idx3 = matrix.get_element(1, 2).unwrap();

    matrix.elements_arena[idx1].real = 1.0;
    matrix.elements_arena[idx2].real = 2.0;

    assert_eq!(
        matrix.get_size(true),
        2,
        "External size should dynamically scale to 2"
    );
    assert_eq!(
        matrix.get_size(false),
        2,
        "Internal size should dynamically scale to 2"
    );

    matrix.clear();
    assert_eq!(
        matrix.elements_arena[idx1].real, 0.0,
        "p_element should be cleared to 0.0"
    );
    assert_eq!(
        matrix.elements_arena[idx2].real, 0.0,
        "p_element should be cleared to 0.0"
    );
}

/// Computes the node voltages of a 2-Resistor circuit:
///
/// Current Source (1A) from Ground(0) to Node 1.
/// R1 (1 Ohm) between Node 1 and Node 2.
/// R2 (1 Ohm) between Node 2 and Ground(0).
///
/// G-Matrix =
/// [  1.0, -1.0 ]
/// [ -1.0,  2.0 ]
///
/// RHS (Current) =
/// [ 1.0, 0.0 ]
///
/// Expected Solution (Voltage):
/// V1 = 2.0 V
/// V2 = 1.0 V
#[test]
fn test_nodal_admittance_solve() {
    let mut matrix = MatrixFrame::new(2, false);

    // Get the element indices for the graph
    let g11 = matrix.get_element(1, 1).unwrap();
    let g12 = matrix.get_element(1, 2).unwrap();
    let g21 = matrix.get_element(2, 1).unwrap();
    let g22 = matrix.get_element(2, 2).unwrap();

    // Load Admittances
    matrix.elements_arena[g11].real += 1.0; // R1
    matrix.elements_arena[g12].real -= 1.0; // R1
    matrix.elements_arena[g21].real -= 1.0; // R1
    matrix.elements_arena[g22].real += 1.0; // R1

    matrix.elements_arena[g22].real += 1.0; // R2

    // Perform factorization
    let error_code = matrix.order_and_factor(None, 0.001, 1e-13, true);
    assert_eq!(
        error_code, 0,
        "Factorization should succeed with spOKAY (0)"
    );

    // Perform Solve
    let rhs = vec![0.0, 1.0, 0.0]; // Node 0 ignored by sparse indexing mappings. Node 1 is 1.0A.
    let mut solution = vec![0.0, 0.0, 0.0];

    matrix.solve(&rhs, &mut solution);

    assert!(
        (solution[1] - 2.0).abs() < EPSILON,
        "V1 Should be 2.0V, got {}",
        solution[1]
    );
    assert!(
        (solution[2] - 1.0).abs() < EPSILON,
        "V2 Should be 1.0V, got {}",
        solution[2]
    );
}
