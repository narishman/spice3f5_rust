use crate::circuit::Circuit;

/// Specifies the numerical integration method
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntegrationMethod {
    Trapezoidal,
    Gear,
}

/// Computes the numerical derivative (current) from charge using the Trapezoidal Rule.
///
/// Returns (geq, ceq) to be stamped into the MNA matrix:
/// * `geq`: Equivalent conductance (dI/dV)
/// * `ceq`: Equivalent current source (-I_history)
/// * `cap`: Current Capacitance (dq/dv)
/// * `state_offset`: The index in the circuit.states array for this device's variable
pub fn integrate(
    circuit: &mut Circuit,
    geq: &mut f64,
    ceq: &mut f64,
    cap: f64,
    state_offset: usize,
) {
    // In SPICE state mapping:
    // state[0] = Current value constructed (e.g. Charge)
    // state[1] = Previous timestep value (e.g. Old Charge)
    // For capacitors, typically offset+0 is Q, offset+1 is I

    // Equivalent of NIintegrate()
    let q0 = circuit.states[0][state_offset]; // Current charge guess
    let q1 = circuit.states[1][state_offset]; // Previous time-step charge
    let i1 = circuit.states[1][state_offset + 1]; // Previous time-step current

    // Trapezoidal rule coefficients: ag[0] = 2/delta, ag[1] = -2/delta
    // circuit.ag array contains standard integration precomputations.
    // However, if we simplify for strict Trapezoidal:
    let two_over_h = 2.0 / circuit.delta;

    *geq = cap * two_over_h;

    let i_history = (two_over_h * q1) + i1; // Expected history impact
    *ceq = -i_history; // The RHS equivalent current source

    // Update current time-step current estimation (so it can be saved after convergence)
    circuit.states[0][state_offset + 1] = (two_over_h * q0) - i_history;
}
