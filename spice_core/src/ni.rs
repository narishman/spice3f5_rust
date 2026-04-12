use crate::circuit::Circuit;
use crate::devices::{Device, DeviceEnum};

pub fn do_dc_op(circuit: &mut Circuit, devices: &mut [DeviceEnum]) -> Result<(), String> {
    // 1. First Pass: Matrix setup and Node bindings
    for dev in devices.iter_mut() {
        if let Err(e) = dev.setup(circuit) {
            return Err(format!("Error setting up device: {}", e));
        }
    }

    // Allocate state vectors based on device setup requirements
    for i in 0..8 {
        circuit.states[i] = vec![0.0; circuit.num_states];
    }

    let mut converged = false;
    let mut iter = 0;

    // 2. Newton-Raphson Iteration Loop
    while iter < circuit.dc_max_iter {
        // Clear matrix and RHS vector
        circuit.matrix.clear();
        circuit.rhs.fill(0.0);
        circuit.noncon = 0;

        // Load conductance and currents into matrix & RHS
        for dev in devices.iter_mut() {
            if let Err(e) = dev.load(circuit) {
                return Err(format!("Error loading device: {}", e));
            }
        }

        // LU Factorization using sparse matrix math
        // In the future: handle matrix singularity exceptions and GMIN stepping here.
        let factor_err = circuit.matrix.order_and_factor(
            Some(&circuit.rhs),
            circuit.reltol,
            circuit.abstol,
            true,
        );
        if factor_err != 0 {
            return Err(format!(
                "Matrix singular or factorization failed. Error code: {}",
                factor_err
            ));
        }

        // Forward and backward substitution to solve for node voltages/currents
        let mut new_states = vec![0.0; circuit.rhs.len()];
        circuit.matrix.solve(&circuit.rhs, &mut new_states);

        // Convergence Test
        // For linear elements (R, V, I), the Newton-Raphson mathematically converges in 1 pass.
        // We evaluate tolerance differences between old states and new.
        let mut all_converged = true;
        if circuit.noncon > 0 {
            all_converged = false;
        }

        for i in 1..new_states.len() {
            let diff = (new_states[i] - circuit.rhs_old[i]).abs();
            let tol = circuit.reltol * new_states[i].abs().max(circuit.rhs_old[i].abs())
                + circuit.volt_tol;

            // If we are iterating non-linear parts, diff > tol fails convergence.
            // Since we implement only basic linear elements now, we simply break early when converged.
            if diff > tol && iter > 0 {
                eprintln!(
                    "Iter {}: Node {} failed convergence. Diff {}, Tol {}, Old {}, New {}",
                    iter, i, diff, tol, circuit.rhs_old[i], new_states[i]
                );
                all_converged = false;
            }
        }

        eprintln!("Iter {}: States evaluated: {:?}", iter, new_states);
        circuit.rhs_old = new_states.clone(); // Correctly track previous solution
        circuit.rhs = new_states.clone(); // Commit computed result to state buffer

        if all_converged && iter > 0 {
            converged = true;
            break;
        }

        // First iteration is always pushed through to derive a starting operating curve.
        iter += 1;
    }

    if !converged {
        return Err(format!(
            "DC Operating Point failed to converge after {} iterations",
            circuit.dc_max_iter
        ));
    }

    Ok(())
}

pub fn do_tran(
    circuit: &mut Circuit,
    devices: &mut [DeviceEnum],
    step_size: f64,
    t_stop: f64,
) -> Result<(), String> {
    // 1. Initial DC Operating Point
    circuit.delta = 0.0;
    do_dc_op(circuit, devices)?;

    // Copy computed DC operating point to RHS_old to start transients
    circuit.rhs_old = circuit.rhs.clone();

    // Copy state 0 into state 1 (starting history)
    circuit.states[1] = circuit.states[0].clone();

    circuit.delta = step_size;
    let mut current_time = 0.0;

    while current_time < t_stop {
        current_time += step_size;

        // Newton-Raphson Loop for this specific time-step
        let mut converged = false;
        let mut iter = 0;

        while iter < circuit.tran_max_iter {
            circuit.matrix.clear();
            circuit.rhs.fill(0.0);
            circuit.noncon = 0;

            for dev in devices.iter_mut() {
                dev.load(circuit)?;
            }

            let factor_err = circuit.matrix.order_and_factor(
                Some(&circuit.rhs),
                circuit.reltol,
                circuit.abstol,
                true,
            );
            if factor_err != 0 {
                return Err(format!("Matrix factorization failed at t={}", current_time));
            }

            let mut new_states = vec![0.0; circuit.rhs.len()];
            circuit.matrix.solve(&circuit.rhs, &mut new_states);

            let mut all_converged = true;
            if circuit.noncon > 0 {
                all_converged = false;
            }

            for i in 1..new_states.len() {
                let diff = (new_states[i] - circuit.rhs_old[i]).abs();
                let tol = circuit.reltol * new_states[i].abs().max(circuit.rhs_old[i].abs())
                    + circuit.volt_tol;
                if diff > tol {
                    all_converged = false;
                }
            }

            circuit.rhs_old = new_states.clone(); // The previous guess for the next NR iteration
            circuit.rhs = new_states.clone();

            if all_converged && iter > 0 {
                converged = true;
                break;
            }
            iter += 1;
        }

        if !converged {
            return Err(format!(
                "Transient failed to converge at t={}",
                current_time
            ));
        }

        // Accept time step: shift states backwards
        circuit.rhs_old = circuit.rhs.clone();
        circuit.states[1] = circuit.states[0].clone();
    }

    Ok(())
}
