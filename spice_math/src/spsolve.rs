//! MATRIX SOLVE MODULE
//!
//! Author: Kenneth S. Kundert (Original C version), Migrated to Rust
//! Advising professor: Alberto Sangiovanni-Vincentelli, UC Berkeley
//!
//! This module contains the forward and backward substitution routines for
//! the sparse matrix solving engine.

use crate::spdefs::MatrixFrame;

impl MatrixFrame {
    /// SOLVE MATRIX EQUATION
    ///
    /// Performs forward elimination and back substitution to find the
    /// unknown vector from the RHS vector and factored matrix. This
    /// routine assumes that the pivots are associated with the lower
    /// triangular (L) matrix and that the diagonal of the upper triangular
    /// (U) matrix consists of ones.
    ///
    /// `rhs` is the input right-hand side.
    /// `solution` is the output data array.
    pub fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) {
        assert!(
            self.factored,
            "Matrix must be factored before calling solve."
        );

        if self.complex {
            self.solve_complex_matrix(rhs, solution);
            return;
        }

        let size = self.size;

        // Scramble / Load intermediate
        // The IntToExtRowMap maps internal row indices back to external ones provided by the user.
        for i in 1..=size {
            let ext_row = self.int_to_ext_row_map[i];
            if ext_row < rhs.len() {
                self.intermediate[i] = rhs[ext_row];
            } else {
                self.intermediate[i] = 0.0;
            }
        }

        // Forward elimination. Solves Lc = b.
        for i in 1..=size {
            let mut temp = self.intermediate[i];

            // This step is skipped if temp equals zero, leveraging sparsity
            if temp != 0.0 {
                let pivot_idx = self.diag[i].expect("Diagonal missing during solve");
                let pivot_real = self.elements_arena[pivot_idx].real;

                temp *= pivot_real;
                self.intermediate[i] = temp;

                let mut p_element = self.elements_arena[pivot_idx].next_in_col;
                while let Some(idx) = p_element {
                    let elem = &self.elements_arena[idx];
                    self.intermediate[elem.row] -= temp * elem.real;
                    p_element = elem.next_in_col;
                }
            }
        }

        // Backward Substitution. Solves Ux = c.
        for i in (1..=size).rev() {
            let mut temp = self.intermediate[i];
            let pivot_idx = self.diag[i].unwrap();

            let mut p_element = self.elements_arena[pivot_idx].next_in_row;
            while let Some(idx) = p_element {
                let elem = &self.elements_arena[idx];
                temp -= elem.real * self.intermediate[elem.col];
                p_element = elem.next_in_row;
            }
            self.intermediate[i] = temp;
        }

        // Unscramble Intermediate vector while placing data into Solution vector.
        for i in 1..=size {
            let ext_col = self.int_to_ext_col_map[i];
            if ext_col < solution.len() {
                solution[ext_col] = self.intermediate[i];
            }
        }
    }

    /// SOLVE TRANSPOSED MATRIX EQUATION
    ///
    /// Performs forward elimination and back substitution to find the
    /// unknown vector from the RHS vector and transposed factored matrix.
    /// This routine is useful when performing sensitivity analysis on a circuit
    /// using the adjoint method.
    pub fn solve_transposed(&mut self, rhs: &[f64], solution: &mut [f64]) {
        assert!(
            self.factored,
            "Matrix must be factored before calling solve_transposed."
        );

        if self.complex {
            self.solve_complex_transposed_matrix(rhs, solution);
            return;
        }

        let size = self.size;

        // Scramble / Load intermediate
        for i in 1..=size {
            let ext_col = self.int_to_ext_col_map[i];
            if ext_col < rhs.len() {
                self.intermediate[i] = rhs[ext_col];
            } else {
                self.intermediate[i] = 0.0;
            }
        }

        // Forward elimination (Transposed)
        for i in 1..=size {
            let temp = self.intermediate[i];
            if temp != 0.0 {
                let pivot_idx = self.diag[i].unwrap();
                let mut p_element = self.elements_arena[pivot_idx].next_in_row;
                while let Some(idx) = p_element {
                    let elem = &self.elements_arena[idx];
                    self.intermediate[elem.col] -= temp * elem.real;
                    p_element = elem.next_in_row;
                }
            }
        }

        // Backward Substitution (Transposed)
        for i in (1..=size).rev() {
            let pivot_idx = self.diag[i].unwrap();
            let mut temp = self.intermediate[i];

            let mut p_element = self.elements_arena[pivot_idx].next_in_col;
            while let Some(idx) = p_element {
                let elem = &self.elements_arena[idx];
                temp -= elem.real * self.intermediate[elem.row];
                p_element = elem.next_in_col;
            }

            self.intermediate[i] = temp * self.elements_arena[pivot_idx].real;
        }

        // Unscramble array into solution
        for i in 1..=size {
            let ext_row = self.int_to_ext_row_map[i];
            if ext_row < solution.len() {
                solution[ext_row] = self.intermediate[i];
            }
        }
    }

    /// SOLVE COMPLEX MATRIX EQUATION (Separated Real/Imaginary vectors, 1:1 with spSEPARATED_COMPLEX_VECTORS)
    pub fn solve_complex(
        &mut self,
        rhs_real: &[f64],
        rhs_imag: &[f64],
        sol_real: &mut [f64],
        sol_imag: &mut [f64],
    ) {
        assert!(
            self.factored,
            "Matrix must be factored before calling solve_complex."
        );

        let size = self.size;
        if self.intermediate.len() < 2 * (size + 1) {
            self.intermediate.resize(2 * (size + 1), 0.0);
        }

        // Initialize Intermediate vector
        for i in 1..=size {
            let ext_row = self.int_to_ext_row_map[i];
            self.intermediate[2 * i] = if ext_row < rhs_real.len() {
                rhs_real[ext_row]
            } else {
                0.0
            };
            self.intermediate[2 * i + 1] = if ext_row < rhs_imag.len() {
                rhs_imag[ext_row]
            } else {
                0.0
            };
        }

        // Forward substitution: Solves Lc = b
        for i in 1..=size {
            let tr = self.intermediate[2 * i];
            let ti = self.intermediate[2 * i + 1];

            if tr != 0.0 || ti != 0.0 {
                let pivot_idx = self.diag[i].expect("Diagonal missing during solve_complex");
                let pr = self.elements_arena[pivot_idx].real;
                let pi = self.elements_arena[pivot_idx].imag;

                // Temp *= Pivot (which stores 1/Pivot)
                let new_tr = tr * pr - ti * pi;
                let new_ti = tr * pi + ti * pr;
                self.intermediate[2 * i] = new_tr;
                self.intermediate[2 * i + 1] = new_ti;

                let mut p_element = self.elements_arena[pivot_idx].next_in_col;
                while let Some(idx) = p_element {
                    let elem = &self.elements_arena[idx];
                    let row = elem.row;
                    self.intermediate[2 * row] -= new_tr * elem.real - new_ti * elem.imag;
                    self.intermediate[2 * row + 1] -= new_tr * elem.imag + new_ti * elem.real;
                    p_element = elem.next_in_col;
                }
            }
        }

        // Backward substitution: Solves Ux = c
        for i in (1..=size).rev() {
            let mut tr = self.intermediate[2 * i];
            let mut ti = self.intermediate[2 * i + 1];
            let pivot_idx = self.diag[i].unwrap();

            let mut p_element = self.elements_arena[pivot_idx].next_in_row;
            while let Some(idx) = p_element {
                let elem = &self.elements_arena[idx];
                let col = elem.col;
                let cr = self.intermediate[2 * col];
                let ci = self.intermediate[2 * col + 1];
                tr -= elem.real * cr - elem.imag * ci;
                ti -= elem.real * ci + elem.imag * cr;
                p_element = elem.next_in_row;
            }
            self.intermediate[2 * i] = tr;
            self.intermediate[2 * i + 1] = ti;
        }

        // Unscramble Intermediate vector while placing data into Solution vectors
        for i in 1..=size {
            let ext_col = self.int_to_ext_col_map[i];
            if ext_col < sol_real.len() {
                sol_real[ext_col] = self.intermediate[2 * i];
            }
            if ext_col < sol_imag.len() {
                sol_imag[ext_col] = self.intermediate[2 * i + 1];
            }
        }
    }

    /// SOLVE COMPLEX TRANSPOSED MATRIX EQUATION (Separated Real/Imaginary vectors)
    pub fn solve_complex_transposed(
        &mut self,
        rhs_real: &[f64],
        rhs_imag: &[f64],
        sol_real: &mut [f64],
        sol_imag: &mut [f64],
    ) {
        assert!(
            self.factored,
            "Matrix must be factored before calling solve_complex_transposed."
        );

        let size = self.size;
        if self.intermediate.len() < 2 * (size + 1) {
            self.intermediate.resize(2 * (size + 1), 0.0);
        }

        for i in 1..=size {
            let ext_col = self.int_to_ext_col_map[i];
            self.intermediate[2 * i] = if ext_col < rhs_real.len() {
                rhs_real[ext_col]
            } else {
                0.0
            };
            self.intermediate[2 * i + 1] = if ext_col < rhs_imag.len() {
                rhs_imag[ext_col]
            } else {
                0.0
            };
        }

        // Forward elimination (Transposed)
        for i in 1..=size {
            let tr = self.intermediate[2 * i];
            let ti = self.intermediate[2 * i + 1];
            if tr != 0.0 || ti != 0.0 {
                let pivot_idx = self.diag[i].unwrap();
                let mut p_element = self.elements_arena[pivot_idx].next_in_row;
                while let Some(idx) = p_element {
                    let elem = &self.elements_arena[idx];
                    let col = elem.col;
                    self.intermediate[2 * col] -= tr * elem.real - ti * elem.imag;
                    self.intermediate[2 * col + 1] -= tr * elem.imag + ti * elem.real;
                    p_element = elem.next_in_row;
                }
            }
        }

        // Backward Substitution (Transposed)
        for i in (1..=size).rev() {
            let pivot_idx = self.diag[i].unwrap();
            let mut tr = self.intermediate[2 * i];
            let mut ti = self.intermediate[2 * i + 1];

            let mut p_element = self.elements_arena[pivot_idx].next_in_col;
            while let Some(idx) = p_element {
                let elem = &self.elements_arena[idx];
                let row = elem.row;
                let rr = self.intermediate[2 * row];
                let ri = self.intermediate[2 * row + 1];
                tr -= elem.real * rr - elem.imag * ri;
                ti -= elem.real * ri + elem.imag * rr;
                p_element = elem.next_in_col;
            }

            let pr = self.elements_arena[pivot_idx].real;
            let pi = self.elements_arena[pivot_idx].imag;
            self.intermediate[2 * i] = tr * pr - ti * pi;
            self.intermediate[2 * i + 1] = tr * pi + ti * pr;
        }

        for i in 1..=size {
            let ext_row = self.int_to_ext_row_map[i];
            if ext_row < sol_real.len() {
                sol_real[ext_row] = self.intermediate[2 * i];
            }
            if ext_row < sol_imag.len() {
                sol_imag[ext_row] = self.intermediate[2 * i + 1];
            }
        }
    }

    /// SOLVE COMPLEX MATRIX EQUATION (Interleaved or real-only fallback)
    fn solve_complex_matrix(&mut self, rhs: &[f64], solution: &mut [f64]) {
        let zeros = vec![0.0; rhs.len()];
        let mut imag_sol = vec![0.0; solution.len()];
        self.solve_complex(rhs, &zeros, solution, &mut imag_sol);
    }

    /// SOLVE COMPLEX TRANSPOSED MATRIX EQUATION (Interleaved or real-only fallback)
    fn solve_complex_transposed_matrix(&mut self, rhs: &[f64], solution: &mut [f64]) {
        let zeros = vec![0.0; rhs.len()];
        let mut imag_sol = vec![0.0; solution.len()];
        self.solve_complex_transposed(rhs, &zeros, solution, &mut imag_sol);
    }
}
