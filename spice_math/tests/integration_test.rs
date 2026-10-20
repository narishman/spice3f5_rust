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
    let _idx3 = matrix.get_element(1, 2).unwrap();

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

/// Tests LU factorization when fill-ins are generated during elimination,
/// as well as re-factorization using `factor()`.
#[test]
fn test_matrix_fillin_and_refactor() {
    // A 4x4 matrix where eliminating row/col 1 creates a fill-in at (2, 3) and (3, 2):
    // [ 4 -1 -1  0 ] [ x1 ]   [ 2 ]
    // [-1  4  0 -1 ] [ x2 ] = [ 2 ]
    // [-1  0  4 -1 ] [ x3 ]   [ 2 ]
    // [ 0 -1 -1  4 ] [ x4 ]   [ 2 ]
    // By symmetry, x1 = x2 = x3 = x4 = 1.0.
    let mut matrix = MatrixFrame::new(4, false);

    let e = |m: &mut MatrixFrame, r: usize, c: usize, val: f64| {
        let idx = m.get_element(r, c).unwrap();
        m.elements_arena[idx].real = val;
    };

    e(&mut matrix, 1, 1, 4.0);
    e(&mut matrix, 1, 2, -1.0);
    e(&mut matrix, 1, 3, -1.0);

    e(&mut matrix, 2, 1, -1.0);
    e(&mut matrix, 2, 2, 4.0);
    e(&mut matrix, 2, 4, -1.0);

    e(&mut matrix, 3, 1, -1.0);
    e(&mut matrix, 3, 3, 4.0);
    e(&mut matrix, 3, 4, -1.0);

    e(&mut matrix, 4, 2, -1.0);
    e(&mut matrix, 4, 3, -1.0);
    e(&mut matrix, 4, 4, 4.0);

    let rhs = vec![0.0, 2.0, 2.0, 2.0, 2.0];
    let err = matrix.order_and_factor(Some(&rhs), 1e-3, 1e-12, true);
    assert_eq!(err, 0);

    let mut sol = vec![0.0; 5];
    matrix.solve(&rhs, &mut sol);
    for i in 1..=4 {
        assert!(
            (sol[i] - 1.0).abs() < EPSILON,
            "x{} should be 1.0, got {}",
            i,
            sol[i]
        );
    }

    // Now test clear + reload + fast `factor()` without reordering
    matrix.clear();
    e(&mut matrix, 1, 1, 4.0);
    e(&mut matrix, 1, 2, -1.0);
    e(&mut matrix, 1, 3, -1.0);
    e(&mut matrix, 2, 1, -1.0);
    e(&mut matrix, 2, 2, 4.0);
    e(&mut matrix, 2, 4, -1.0);
    e(&mut matrix, 3, 1, -1.0);
    e(&mut matrix, 3, 3, 4.0);
    e(&mut matrix, 3, 4, -1.0);
    e(&mut matrix, 4, 2, -1.0);
    e(&mut matrix, 4, 3, -1.0);
    e(&mut matrix, 4, 4, 4.0);

    matrix.partition(PartitionMode::AutoPartition);
    let err2 = matrix.factor();
    assert_eq!(err2, 0);

    let rhs2 = vec![0.0, 4.0, 4.0, 4.0, 4.0];
    let mut sol2 = vec![0.0; 5];
    matrix.solve(&rhs2, &mut sol2);
    for i in 1..=4 {
        assert!(
            (sol2[i] - 2.0).abs() < EPSILON,
            "x{} should be 2.0, got {}",
            i,
            sol2[i]
        );
    }
}

/// Tests complex matrix LU factorization and solving.
#[test]
fn test_complex_matrix_solve() {
    // 2x2 complex system:
    // (1 + j1) x1 + (-1) x2 = 1 + j0
    // (-1) x1 + (2 + j1) x2 = 0 + j0
    // Det = (1+j)(2+j) - 1 = (2 + 3j - 1) - 1 = 3j
    // x1 = (2+j) / (3j) = (1 - 2j) / 3 = 1/3 - j 2/3
    // x2 = 1 / (3j) = -j / 3 = 0 - j 1/3
    let mut matrix = MatrixFrame::new(2, true);

    let g11 = matrix.get_element(1, 1).unwrap();
    let g12 = matrix.get_element(1, 2).unwrap();
    let g21 = matrix.get_element(2, 1).unwrap();
    let g22 = matrix.get_element(2, 2).unwrap();

    matrix.elements_arena[g11].real = 1.0;
    matrix.elements_arena[g11].imag = 1.0;
    matrix.elements_arena[g12].real = -1.0;
    matrix.elements_arena[g21].real = -1.0;
    matrix.elements_arena[g22].real = 2.0;
    matrix.elements_arena[g22].imag = 1.0;

    let err = matrix.order_and_factor(None, 1e-3, 1e-12, true);
    assert_eq!(err, 0);

    let rhs_r = vec![0.0, 1.0, 0.0];
    let rhs_i = vec![0.0, 0.0, 0.0];
    let mut sol_r = vec![0.0; 3];
    let mut sol_i = vec![0.0; 3];

    matrix.solve_complex(&rhs_r, &rhs_i, &mut sol_r, &mut sol_i);

    assert!((sol_r[1] - 1.0 / 3.0).abs() < EPSILON);
    assert!((sol_i[1] - (-2.0 / 3.0)).abs() < EPSILON);
    assert!((sol_r[2] - 0.0).abs() < EPSILON);
    assert!((sol_i[2] - (-1.0 / 3.0)).abs() < EPSILON);
}
