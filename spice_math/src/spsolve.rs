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

    /// SOLVE COMPLEX MATRIX EQUATION
    fn solve_complex_matrix(&mut self, rhs: &[f64], solution: &mut [f64]) {
        // Translation stub for complex factorization output evaluation
    }

    /// SOLVE COMPLEX TRANSPOSED MATRIX EQUATION
    fn solve_complex_transposed_matrix(&mut self, rhs: &[f64], solution: &mut [f64]) {
        // Translation stub for complex transposed factorization output evaluation
    }
}
