//! MATRIX FACTORIZATION MODULE
//!
//! Author: Kenneth S. Kundert (Original C version), Migrated to Rust
//! Advising professor: Alberto Sangiovanni-Vincentelli, UC Berkeley
//!
//! This module contains the routines to factor the matrix into LU form using
//! the Markowitz pivot strategy customized for sparse circuit matrices.

use crate::spdefs::MatrixFrame;

#[derive(Debug, PartialEq, Eq)]
pub enum PartitionMode {
    DefaultPartition,
    DirectPartition,
    IndirectPartition,
    AutoPartition,
}

impl MatrixFrame {
    /// ORDER AND FACTOR MATRIX
    ///
    /// This routine chooses a pivot order for the matrix and factors it
    /// into LU form using Gauss's method where the pivots are associated
    /// with L and the diagonal terms of U are one.
    pub fn order_and_factor(
        &mut self,
        rhs: Option<&[f64]>,
        rel_threshold: f64,
        abs_threshold: f64,
        diag_pivoting: bool,
    ) -> i32 {
        if self.factored {
            return self.error; // Should not be factored before calling
        }

        self.error = 0; // spOKAY

        let mut rel_t = rel_threshold;
        if rel_t <= 0.0 || rel_t > 1.0 {
            rel_t = self.rel_threshold;
        }
        self.rel_threshold = rel_t;

        let mut abs_t = abs_threshold;
        if abs_t < 0.0 {
            abs_t = self.abs_threshold;
        }
        self.abs_threshold = abs_t;

        let mut reordering_required = false;
        let mut step = 1;

        if !self.needs_ordering {
            // Already factored once, re-ordering not needed
            for s in 1..=self.size {
                let p_pivot = self.diag[s];
                if p_pivot.is_none() {
                    reordering_required = true;
                    break;
                }

                let pivot_idx = p_pivot.unwrap();
                let pivot_mag = self.element_mag(pivot_idx);
                let largest_in_col =
                    self.find_largest_in_col(self.elements_arena[pivot_idx].next_in_col);

                if largest_in_col * self.rel_threshold < pivot_mag {
                    if self.complex {
                        self.complex_row_col_elimination(pivot_idx);
                    } else {
                        self.real_row_col_elimination(pivot_idx);
                    }
                } else {
                    reordering_required = true;
                    step = s;
                    break;
                }
            }

            if !reordering_required {
                self.needs_ordering = false;
                self.reordered = true;
                self.factored = true;
                return self.error;
            }
        } else {
            // First time factorization (or matrix structure changed)
            if !self.rows_linked {
                self.link_rows();
            }
            self.create_internal_vectors();
            if self.error >= 2 {
                // spFATAL
                return self.error;
            }
        }

        // Form initial Markowitz products
        self.count_markowitz(rhs, step);
        self.markowitz_products(step);
        self.max_row_count_in_lower_tri = 0;

        for s in step..=self.size {
            let p_pivot = self.search_for_pivot(s, diag_pivoting);

            if p_pivot.is_none() {
                return self.matrix_is_singular(s);
            }

            let pivot_idx = p_pivot.unwrap();
            self.exchange_rows_and_cols(pivot_idx, s);

            if self.complex {
                self.complex_row_col_elimination(pivot_idx);
            } else {
                self.real_row_col_elimination(pivot_idx);
            }

            if self.error >= 2 {
                return self.error;
            }

            self.update_markowitz_numbers(pivot_idx);
        }

        self.needs_ordering = false;
        self.reordered = true;
        self.factored = true;

        self.error
    }

    /// FACTOR MATRIX
    ///
    /// Companion routine to `order_and_factor`. Cannot change ordering.
    pub fn factor(&mut self) -> i32 {
        if self.factored {
            return self.error;
        }

        if self.needs_ordering {
            return self.order_and_factor(None, 0.0, 0.0, true);
        }

        if !self.partitioned {
            self.partition(PartitionMode::DefaultPartition);
        }

        if self.complex {
            return self.factor_complex_matrix();
        }

        let pivot_idx = match self.diag[1] {
            Some(idx) => idx,
            None => return self.zero_pivot(1),
        };
        if self.elements_arena[pivot_idx].real == 0.0 {
            return self.zero_pivot(1);
        }
        self.elements_arena[pivot_idx].real = 1.0 / self.elements_arena[pivot_idx].real;

        for step in 2..=self.size {
            // In safe Rust using contiguous arena vectors, both direct and indirect
            // partition modes execute via the intermediate scatter-update-gather buffer.
            let mut p_elem_idx = self.first_in_col[step];
            while let Some(idx) = p_elem_idx {
                let row = self.elements_arena[idx].row;
                self.intermediate[row] = self.elements_arena[idx].real;
                p_elem_idx = self.elements_arena[idx].next_in_col;
            }

            // Update
            let mut p_col_idx = self.first_in_col[step];
            while let Some(col_idx) = p_col_idx {
                let col_row = self.elements_arena[col_idx].row;
                if col_row >= step {
                    break;
                }

                let p_diag_idx = self.diag[col_row].unwrap();
                let diag_real = self.elements_arena[p_diag_idx].real;

                self.elements_arena[col_idx].real = self.intermediate[col_row] * diag_real;
                let mult = self.elements_arena[col_idx].real;

                let mut p_elem2_idx = self.elements_arena[p_diag_idx].next_in_col;
                while let Some(idx2) = p_elem2_idx {
                    let row2 = self.elements_arena[idx2].row;
                    self.intermediate[row2] -= mult * self.elements_arena[idx2].real;
                    p_elem2_idx = self.elements_arena[idx2].next_in_col;
                }

                p_col_idx = self.elements_arena[col_idx].next_in_col;
            }

            // Gather
            let step_diag_idx = match self.diag[step] {
                Some(idx) => idx,
                None => return self.zero_pivot(step),
            };
            let mut p_elem3_idx = self.elements_arena[step_diag_idx].next_in_col;
            while let Some(idx3) = p_elem3_idx {
                let row3 = self.elements_arena[idx3].row;
                self.elements_arena[idx3].real = self.intermediate[row3];
                p_elem3_idx = self.elements_arena[idx3].next_in_col;
            }

            if self.intermediate[step] == 0.0 {
                return self.zero_pivot(step);
            }

            self.elements_arena[step_diag_idx].real = 1.0 / self.intermediate[step];
        }

        self.factored = true;
        self.error = 0;
        self.error
    }

    /// FACTOR COMPLEX MATRIX
    fn factor_complex_matrix(&mut self) -> i32 {
        let pivot_idx = match self.diag[1] {
            Some(idx) => idx,
            None => return self.zero_pivot(1),
        };
        if self.element_mag(pivot_idx) == 0.0 {
            return self.zero_pivot(1);
        }
        let (pr, pi) = Self::cmplx_reciprocal(
            self.elements_arena[pivot_idx].real,
            self.elements_arena[pivot_idx].imag,
        );
        self.elements_arena[pivot_idx].real = pr;
        self.elements_arena[pivot_idx].imag = pi;

        // Ensure intermediate has space for 2 * (size + 1)
        if self.intermediate.len() < 2 * (self.size + 1) {
            self.intermediate.resize(2 * (self.size + 1), 0.0);
        }

        for step in 2..=self.size {
            // Scatter complex values into self.intermediate[2*row] and self.intermediate[2*row + 1]
            let mut p_elem_idx = self.first_in_col[step];
            while let Some(idx) = p_elem_idx {
                let row = self.elements_arena[idx].row;
                self.intermediate[2 * row] = self.elements_arena[idx].real;
                self.intermediate[2 * row + 1] = self.elements_arena[idx].imag;
                p_elem_idx = self.elements_arena[idx].next_in_col;
            }

            // Update column
            let mut p_col_idx = self.first_in_col[step];
            while let Some(col_idx) = p_col_idx {
                let col_row = self.elements_arena[col_idx].row;
                if col_row >= step {
                    break;
                }

                let p_diag_idx = self.diag[col_row].unwrap();
                let dr = self.elements_arena[p_diag_idx].real;
                let di = self.elements_arena[p_diag_idx].imag;
                let sr = self.intermediate[2 * col_row];
                let si = self.intermediate[2 * col_row + 1];

                // Mult = Dest[col_row] * Diag[col_row]
                let mult_r = sr * dr - si * di;
                let mult_i = sr * di + si * dr;
                self.elements_arena[col_idx].real = mult_r;
                self.elements_arena[col_idx].imag = mult_i;

                let mut p_elem2_idx = self.elements_arena[p_diag_idx].next_in_col;
                while let Some(idx2) = p_elem2_idx {
                    let row2 = self.elements_arena[idx2].row;
                    let er = self.elements_arena[idx2].real;
                    let ei = self.elements_arena[idx2].imag;
                    self.intermediate[2 * row2] -= mult_r * er - mult_i * ei;
                    self.intermediate[2 * row2 + 1] -= mult_r * ei + mult_i * er;
                    p_elem2_idx = self.elements_arena[idx2].next_in_col;
                }

                p_col_idx = self.elements_arena[col_idx].next_in_col;
            }

            // Gather
            let step_diag_idx = match self.diag[step] {
                Some(idx) => idx,
                None => return self.zero_pivot(step),
            };
            let mut p_elem3_idx = self.elements_arena[step_diag_idx].next_in_col;
            while let Some(idx3) = p_elem3_idx {
                let row3 = self.elements_arena[idx3].row;
                self.elements_arena[idx3].real = self.intermediate[2 * row3];
                self.elements_arena[idx3].imag = self.intermediate[2 * row3 + 1];
                p_elem3_idx = self.elements_arena[idx3].next_in_col;
            }

            let piv_r = self.intermediate[2 * step];
            let piv_i = self.intermediate[2 * step + 1];
            if piv_r.abs() + piv_i.abs() == 0.0 {
                return self.zero_pivot(step);
            }
            let (rec_r, rec_i) = Self::cmplx_reciprocal(piv_r, piv_i);
            self.elements_arena[step_diag_idx].real = rec_r;
            self.elements_arena[step_diag_idx].imag = rec_i;
        }

        self.factored = true;
        self.error = 0;
        self.error
    }

    /// CREATE INTERNAL VECTORS
    fn create_internal_vectors(&mut self) {
        self.markowitz_row.resize(self.size + 1, 0);
        self.markowitz_col.resize(self.size + 1, 0);
        self.markowitz_prod.resize(self.size + 2, 0);

        self.do_real_direct.resize(self.size + 1, false);
        self.do_cmplx_direct.resize(self.size + 1, false);

        let inter_size = if self.complex {
            2 * (self.size + 1)
        } else {
            self.size + 1
        };
        self.intermediate.resize(inter_size, 0.0);
        self.internal_vectors_allocated = true;
    }

    /// PARTITION MATRIX
    pub fn partition(&mut self, mode: PartitionMode) {
        if self.partitioned {
            return;
        }

        self.do_real_direct.resize(self.size + 1, false);
        self.do_cmplx_direct.resize(self.size + 1, false);
        self.partitioned = true;

        let effective_mode = match mode {
            PartitionMode::DefaultPartition => PartitionMode::AutoPartition,
            other => other,
        };

        match effective_mode {
            PartitionMode::DirectPartition => {
                for step in 1..=self.size {
                    self.do_real_direct[step] = true;
                    self.do_cmplx_direct[step] = true;
                }
                return;
            }
            PartitionMode::IndirectPartition => {
                for step in 1..=self.size {
                    self.do_real_direct[step] = false;
                    self.do_cmplx_direct[step] = false;
                }
                return;
            }
            PartitionMode::AutoPartition | PartitionMode::DefaultPartition => {}
        }

        // Count all operations needed in factoring matrix (1:1 with spPartition)
        for step in 1..=self.size {
            let mut nc = 0usize;
            let mut no = 0usize;
            let mut nm = 0usize;

            let mut p_element = self.first_in_col[step];
            while let Some(idx) = p_element {
                nc += 1;
                p_element = self.elements_arena[idx].next_in_col;
            }

            let mut p_column = self.first_in_col[step];
            while let Some(col_idx) = p_column {
                let col_row = self.elements_arena[col_idx].row;
                if col_row >= step {
                    break;
                }
                if let Some(diag_idx) = self.diag[col_row] {
                    nm += 1;
                    let mut p_sub = self.elements_arena[diag_idx].next_in_col;
                    while let Some(sub_idx) = p_sub {
                        no += 1;
                        p_sub = self.elements_arena[sub_idx].next_in_col;
                    }
                }
                p_column = self.elements_arena[col_idx].next_in_col;
            }

            self.do_real_direct[step] = (nm + no) as isize > (3 * nc as isize - 2 * nm as isize);
            self.do_cmplx_direct[step] = (nm + no) as isize > (7 * nc as isize - 4 * nm as isize);
        }
    }

    fn search_entire_matrix(&mut self, step: usize) -> Option<usize> {
        let mut chosen_pivot = None;
        let mut largest_element_mag = 0.0;
        let mut p_largest_element = None;
        let mut min_markowitz_product = i64::MAX;
        let mut number_of_ties = 0;
        let mut ratio_of_accepted = 0.0;

        for i in step..=self.size {
            let mut p_element = self.first_in_col[i];

            while let Some(idx) = p_element {
                if self.elements_arena[idx].row >= step {
                    break;
                }
                p_element = self.elements_arena[idx].next_in_col;
            }

            let largest_in_col = self.find_largest_in_col(p_element);
            if largest_in_col == 0.0 {
                continue;
            }

            while let Some(idx) = p_element {
                let magnitude = self.element_mag(idx);
                if magnitude > largest_element_mag {
                    largest_element_mag = magnitude;
                    p_largest_element = Some(idx);
                }

                let row = self.elements_arena[idx].row;
                let col = self.elements_arena[idx].col;
                let product = (self.markowitz_row[row] as i64) * (self.markowitz_col[col] as i64);

                if product <= min_markowitz_product
                    && magnitude > self.rel_threshold * largest_in_col
                    && magnitude > self.abs_threshold
                {
                    if product < min_markowitz_product {
                        chosen_pivot = Some(idx);
                        min_markowitz_product = product;
                        ratio_of_accepted = largest_in_col / magnitude;
                        number_of_ties = 0;
                    } else {
                        number_of_ties += 1;
                        let ratio = largest_in_col / magnitude;
                        if ratio < ratio_of_accepted {
                            chosen_pivot = Some(idx);
                            ratio_of_accepted = ratio;
                        }
                        // TIES_MULTIPLIER = 5
                        if number_of_ties as i64 >= min_markowitz_product * 5 {
                            return chosen_pivot;
                        }
                    }
                }
                p_element = self.elements_arena[idx].next_in_col;
            }
        }

        if chosen_pivot.is_some() {
            return chosen_pivot;
        }

        if largest_element_mag == 0.0 {
            self.error = 1; // spSINGULAR
            return None;
        }

        self.error = 3; // spSMALL_PIVOT
        p_largest_element
    }

    fn matrix_is_singular(&mut self, step: usize) -> i32 {
        self.error = 1; // spSINGULAR
        self.singular_row = self.int_to_ext_row_map[step];
        self.singular_col = self.int_to_ext_col_map[step];
        1
    }

    fn exchange_rows_and_cols(&mut self, pivot_idx: usize, step: usize) {
        let p_pivot = pivot_idx;
        let row = self.elements_arena[p_pivot].row;
        let col = self.elements_arena[p_pivot].col;
        self.pivots_original_row = row;
        self.pivots_original_col = col;

        if row == step && col == step {
            return;
        }

        if row == col {
            self.exchange_rows(step, row);
            self.exchange_cols(step, col);
            self.markowitz_prod.swap(step, row);
            self.diag.swap(row, step);
        } else {
            let old_prod_step = self.markowitz_prod[step];
            let old_prod_row = self.markowitz_prod[row];
            let old_prod_col = self.markowitz_prod[col];

            if row != step {
                self.exchange_rows(step, row);
                self.number_of_interchanges_is_odd = !self.number_of_interchanges_is_odd;
                self.markowitz_prod[row] = self.markowitz_row[row] * self.markowitz_col[row];
                if (self.markowitz_prod[row] == 0) != (old_prod_row == 0) {
                    if old_prod_row == 0 {
                        self.singletons -= 1;
                    } else {
                        self.singletons += 1;
                    }
                }
            }

            if col != step {
                self.exchange_cols(step, col);
                self.number_of_interchanges_is_odd = !self.number_of_interchanges_is_odd;
                self.markowitz_prod[col] = self.markowitz_row[col] * self.markowitz_col[col];
                if (self.markowitz_prod[col] == 0) != (old_prod_col == 0) {
                    if old_prod_col == 0 {
                        self.singletons -= 1;
                    } else {
                        self.singletons += 1;
                    }
                }
                self.diag[col] = self.find_element_in_col(col, col, false);
            }

            if row != step {
                self.diag[row] = self.find_element_in_col(row, row, false);
            }
            self.diag[step] = self.find_element_in_col(step, step, false);

            self.markowitz_prod[step] = self.markowitz_row[step] * self.markowitz_col[step];
            if (self.markowitz_prod[step] == 0) != (old_prod_step == 0) {
                if old_prod_step == 0 {
                    self.singletons -= 1;
                } else {
                    self.singletons += 1;
                }
            }
        }
    }

    fn exchange_rows(&mut self, mut row1: usize, mut row2: usize) {
        if row1 > row2 {
            std::mem::swap(&mut row1, &mut row2);
        }
        let mut p_row1 = self.first_in_row[row1];
        let mut p_row2 = self.first_in_row[row2];

        while p_row1.is_some() || p_row2.is_some() {
            let col;
            let mut e1 = None;
            let mut e2 = None;

            if p_row1.is_none() {
                let r2 = p_row2.unwrap();
                col = self.elements_arena[r2].col;
                e2 = p_row2;
                p_row2 = self.elements_arena[r2].next_in_row;
            } else if p_row2.is_none() {
                let r1 = p_row1.unwrap();
                col = self.elements_arena[r1].col;
                e1 = p_row1;
                p_row1 = self.elements_arena[r1].next_in_row;
            } else {
                let r1 = p_row1.unwrap();
                let r2 = p_row2.unwrap();
                let col1 = self.elements_arena[r1].col;
                let col2 = self.elements_arena[r2].col;

                if col1 < col2 {
                    col = col1;
                    e1 = p_row1;
                    p_row1 = self.elements_arena[r1].next_in_row;
                } else if col1 > col2 {
                    col = col2;
                    e2 = p_row2;
                    p_row2 = self.elements_arena[r2].next_in_row;
                } else {
                    col = col1;
                    e1 = p_row1;
                    e2 = p_row2;
                    p_row1 = self.elements_arena[r1].next_in_row;
                    p_row2 = self.elements_arena[r2].next_in_row;
                }
            }
            self.exchange_col_elements(row1, e1, row2, e2, col);
        }

        self.markowitz_row.swap(row1, row2);
        self.first_in_row.swap(row1, row2);
        self.int_to_ext_row_map.swap(row1, row2);
        let ext1 = self.int_to_ext_row_map[row1];
        let ext2 = self.int_to_ext_row_map[row2];
        self.ext_to_int_row_map[ext1] = row1;
        self.ext_to_int_row_map[ext2] = row2;
    }

    fn exchange_cols(&mut self, mut col1: usize, mut col2: usize) {
        if col1 > col2 {
            std::mem::swap(&mut col1, &mut col2);
        }
        let mut p_col1 = self.first_in_col[col1];
        let mut p_col2 = self.first_in_col[col2];

        while p_col1.is_some() || p_col2.is_some() {
            let row;
            let mut e1 = None;
            let mut e2 = None;

            if p_col1.is_none() {
                let c2 = p_col2.unwrap();
                row = self.elements_arena[c2].row;
                e2 = p_col2;
                p_col2 = self.elements_arena[c2].next_in_col;
            } else if p_col2.is_none() {
                let c1 = p_col1.unwrap();
                row = self.elements_arena[c1].row;
                e1 = p_col1;
                p_col1 = self.elements_arena[c1].next_in_col;
            } else {
                let c1 = p_col1.unwrap();
                let c2 = p_col2.unwrap();
                let r1 = self.elements_arena[c1].row;
                let r2 = self.elements_arena[c2].row;

                if r1 < r2 {
                    row = r1;
                    e1 = p_col1;
                    p_col1 = self.elements_arena[c1].next_in_col;
                } else if r1 > r2 {
                    row = r2;
                    e2 = p_col2;
                    p_col2 = self.elements_arena[c2].next_in_col;
                } else {
                    row = r1;
                    e1 = p_col1;
                    e2 = p_col2;
                    p_col1 = self.elements_arena[c1].next_in_col;
                    p_col2 = self.elements_arena[c2].next_in_col;
                }
            }
            self.exchange_row_elements(col1, e1, col2, e2, row);
        }

        self.markowitz_col.swap(col1, col2);
        self.first_in_col.swap(col1, col2);
        self.int_to_ext_col_map.swap(col1, col2);
        let ext1 = self.int_to_ext_col_map[col1];
        let ext2 = self.int_to_ext_col_map[col2];
        self.ext_to_int_col_map[ext1] = col1;
        self.ext_to_int_col_map[ext2] = col2;
    }

    fn exchange_col_elements(
        &mut self,
        row1: usize,
        e1: Option<usize>,
        row2: usize,
        e2: Option<usize>,
        col: usize,
    ) {
        let mut prev1 = None;
        let mut curr = self.first_in_col[col];
        while let Some(c) = curr {
            if self.elements_arena[c].row < row1 {
                prev1 = Some(c);
                curr = self.elements_arena[c].next_in_col;
            } else {
                break;
            }
        }

        if let Some(e1_idx) = e1 {
            let below1 = self.elements_arena[e1_idx].next_in_col;
            if e2.is_none() {
                if let Some(b1) = below1 {
                    if self.elements_arena[b1].row < row2 {
                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_col = below1;
                        } else {
                            self.first_in_col[col] = below1;
                        }

                        let mut p2_prev = None;
                        let mut p2_curr = below1;
                        while let Some(c) = p2_curr {
                            if self.elements_arena[c].row < row2 {
                                p2_prev = Some(c);
                                p2_curr = self.elements_arena[c].next_in_col;
                            } else {
                                break;
                            }
                        }

                        if let Some(p2) = p2_prev {
                            self.elements_arena[p2].next_in_col = Some(e1_idx);
                        } else {
                            self.first_in_col[col] = Some(e1_idx);
                        }

                        self.elements_arena[e1_idx].next_in_col = p2_curr;
                    }
                }
                self.elements_arena[e1_idx].row = row2;
            } else {
                let e2_idx = e2.unwrap();
                if let Some(b1) = below1 {
                    if self.elements_arena[b1].row == row2 {
                        self.elements_arena[e1_idx].next_in_col =
                            self.elements_arena[e2_idx].next_in_col;
                        self.elements_arena[e2_idx].next_in_col = Some(e1_idx);
                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_col = Some(e2_idx);
                        } else {
                            self.first_in_col[col] = Some(e2_idx);
                        }
                    } else {
                        let mut p2_prev = None;
                        let mut p2_curr = below1;
                        while let Some(c) = p2_curr {
                            if self.elements_arena[c].row < row2 {
                                p2_prev = Some(c);
                                p2_curr = self.elements_arena[c].next_in_col;
                            } else {
                                break;
                            }
                        }
                        let below2 = self.elements_arena[e2_idx].next_in_col;

                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_col = Some(e2_idx);
                        } else {
                            self.first_in_col[col] = Some(e2_idx);
                        }
                        self.elements_arena[e2_idx].next_in_col = below1;

                        if let Some(p2) = p2_prev {
                            self.elements_arena[p2].next_in_col = Some(e1_idx);
                        } else {
                            self.first_in_col[col] = Some(e1_idx);
                        }
                        self.elements_arena[e1_idx].next_in_col = below2;
                    }
                }
                self.elements_arena[e1_idx].row = row2;
                self.elements_arena[e2_idx].row = row1;
            }
        } else {
            let e2_idx = e2.unwrap();
            let mut p2_prev = None;
            let mut p2_curr = curr;
            while let Some(c) = p2_curr {
                if self.elements_arena[c].row < row2 {
                    p2_prev = Some(c);
                    p2_curr = self.elements_arena[c].next_in_col;
                } else {
                    break;
                }
            }
            if let Some(curr_idx) = curr {
                if self.elements_arena[curr_idx].row != row2 {
                    let below2 = self.elements_arena[e2_idx].next_in_col;
                    if let Some(p2) = p2_prev {
                        self.elements_arena[p2].next_in_col = below2;
                    } else {
                        self.first_in_col[col] = below2;
                    }

                    if let Some(p1) = prev1 {
                        self.elements_arena[p1].next_in_col = Some(e2_idx);
                    } else {
                        self.first_in_col[col] = Some(e2_idx);
                    }
                    self.elements_arena[e2_idx].next_in_col = curr;
                }
            }
            self.elements_arena[e2_idx].row = row1;
        }
    }

    fn exchange_row_elements(
        &mut self,
        col1: usize,
        e1: Option<usize>,
        col2: usize,
        e2: Option<usize>,
        row: usize,
    ) {
        let mut prev1 = None;
        let mut curr = self.first_in_row[row];
        while let Some(c) = curr {
            if self.elements_arena[c].col < col1 {
                prev1 = Some(c);
                curr = self.elements_arena[c].next_in_row;
            } else {
                break;
            }
        }

        if let Some(e1_idx) = e1 {
            let right1 = self.elements_arena[e1_idx].next_in_row;
            if e2.is_none() {
                if let Some(r1) = right1 {
                    if self.elements_arena[r1].col < col2 {
                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_row = right1;
                        } else {
                            self.first_in_row[row] = right1;
                        }

                        let mut p2_prev = None;
                        let mut p2_curr = right1;
                        while let Some(c) = p2_curr {
                            if self.elements_arena[c].col < col2 {
                                p2_prev = Some(c);
                                p2_curr = self.elements_arena[c].next_in_row;
                            } else {
                                break;
                            }
                        }

                        if let Some(p2) = p2_prev {
                            self.elements_arena[p2].next_in_row = Some(e1_idx);
                        } else {
                            self.first_in_row[row] = Some(e1_idx);
                        }
                        self.elements_arena[e1_idx].next_in_row = p2_curr;
                    }
                }
                self.elements_arena[e1_idx].col = col2;
            } else {
                let e2_idx = e2.unwrap();
                if let Some(r1) = right1 {
                    if self.elements_arena[r1].col == col2 {
                        self.elements_arena[e1_idx].next_in_row =
                            self.elements_arena[e2_idx].next_in_row;
                        self.elements_arena[e2_idx].next_in_row = Some(e1_idx);
                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_row = Some(e2_idx);
                        } else {
                            self.first_in_row[row] = Some(e2_idx);
                        }
                    } else {
                        let mut p2_prev = None;
                        let mut p2_curr = right1;
                        while let Some(c) = p2_curr {
                            if self.elements_arena[c].col < col2 {
                                p2_prev = Some(c);
                                p2_curr = self.elements_arena[c].next_in_row;
                            } else {
                                break;
                            }
                        }
                        let right2 = self.elements_arena[e2_idx].next_in_row;

                        if let Some(p1) = prev1 {
                            self.elements_arena[p1].next_in_row = Some(e2_idx);
                        } else {
                            self.first_in_row[row] = Some(e2_idx);
                        }
                        self.elements_arena[e2_idx].next_in_row = right1;

                        if let Some(p2) = p2_prev {
                            self.elements_arena[p2].next_in_row = Some(e1_idx);
                        } else {
                            self.first_in_row[row] = Some(e1_idx);
                        }
                        self.elements_arena[e1_idx].next_in_row = right2;
                    }
                }
                self.elements_arena[e1_idx].col = col2;
                self.elements_arena[e2_idx].col = col1;
            }
        } else {
            let e2_idx = e2.unwrap();
            let mut p2_prev = None;
            let mut p2_curr = curr;
            while let Some(c) = p2_curr {
                if self.elements_arena[c].col < col2 {
                    p2_prev = Some(c);
                    p2_curr = self.elements_arena[c].next_in_row;
                } else {
                    break;
                }
            }
            if let Some(curr_idx) = curr {
                if self.elements_arena[curr_idx].col != col2 {
                    let right2 = self.elements_arena[e2_idx].next_in_row;
                    if let Some(p2) = p2_prev {
                        self.elements_arena[p2].next_in_row = right2;
                    } else {
                        self.first_in_row[row] = right2;
                    }

                    if let Some(p1) = prev1 {
                        self.elements_arena[p1].next_in_row = Some(e2_idx);
                    } else {
                        self.first_in_row[row] = Some(e2_idx);
                    }
                    self.elements_arena[e2_idx].next_in_row = curr;
                }
            }
            self.elements_arena[e2_idx].col = col1;
        }
    }

    fn update_markowitz_numbers(&mut self, pivot_idx: usize) {
        let mut col_ptr = self.elements_arena[pivot_idx].next_in_col;
        while let Some(c_idx) = col_ptr {
            let row = self.elements_arena[c_idx].row;
            self.markowitz_row[row] = self.markowitz_row[row].saturating_sub(1);
            self.markowitz_prod[row] =
                self.markowitz_row[row].saturating_mul(self.markowitz_col[row]);
            if self.markowitz_row[row] == 0 {
                self.singletons += 1;
            }
            col_ptr = self.elements_arena[c_idx].next_in_col;
        }

        let mut row_ptr = self.elements_arena[pivot_idx].next_in_row;
        while let Some(r_idx) = row_ptr {
            let col = self.elements_arena[r_idx].col;
            self.markowitz_col[col] = self.markowitz_col[col].saturating_sub(1);
            self.markowitz_prod[col] =
                self.markowitz_row[col].saturating_mul(self.markowitz_col[col]);
            if self.markowitz_col[col] == 0 && self.markowitz_row[col] != 0 {
                self.singletons += 1;
            }
            row_ptr = self.elements_arena[r_idx].next_in_row;
        }
    }

    fn create_fillin(&mut self, row: usize, col: usize) -> Option<usize> {
        let mut p_element = self.first_in_col[col];
        let mut prev_idx = None;
        while let Some(idx) = p_element {
            if self.elements_arena[idx].row < row {
                prev_idx = Some(idx);
                p_element = self.elements_arena[idx].next_in_col;
            } else {
                break;
            }
        }

        let fillin_idx = self.create_element(row, col, prev_idx, true)?;

        // Update Markowitz counts and products (1:1 with CreateFillin in spfactor.c)
        self.markowitz_row[row] += 1;
        self.markowitz_prod[row] = self.markowitz_row[row].saturating_mul(self.markowitz_col[row]);
        if self.markowitz_row[row] == 1 && self.markowitz_col[row] != 0 {
            self.singletons = self.singletons.saturating_sub(1);
        }

        self.markowitz_col[col] += 1;
        self.markowitz_prod[col] = self.markowitz_col[col].saturating_mul(self.markowitz_row[col]);
        if self.markowitz_row[col] != 0 && self.markowitz_col[col] == 1 {
            self.singletons = self.singletons.saturating_sub(1);
        }

        Some(fillin_idx)
    }

    fn link_rows(&mut self) {
        // Build row links from column links
        for col in 1..=self.size {
            let mut elem_idx = self.first_in_col[col];
            while let Some(idx) = elem_idx {
                let row = self.elements_arena[idx].row;
                // Insert into row
                let mut p_curr_row = self.first_in_row[row];
                let mut p_prev_row = None;

                while let Some(curr_idx) = p_curr_row {
                    if self.elements_arena[curr_idx].col < col {
                        p_prev_row = Some(curr_idx);
                        p_curr_row = self.elements_arena[curr_idx].next_in_row;
                    } else {
                        break;
                    }
                }

                if let Some(prev_idx) = p_prev_row {
                    self.elements_arena[idx].next_in_row =
                        self.elements_arena[prev_idx].next_in_row;
                    self.elements_arena[prev_idx].next_in_row = Some(idx);
                } else {
                    self.elements_arena[idx].next_in_row = self.first_in_row[row];
                    self.first_in_row[row] = Some(idx);
                }

                elem_idx = self.elements_arena[idx].next_in_col;
            }
        }
        self.rows_linked = true;
    }

    fn zero_pivot(&mut self, step: usize) -> i32 {
        self.error = 2; // spZERO_DIAG
        self.singular_row = self.int_to_ext_row_map[step];
        self.singular_col = self.int_to_ext_col_map[step];
        2
    }

    pub(crate) fn cmplx_reciprocal(den_r: f64, den_i: f64) -> (f64, f64) {
        if (den_r >= den_i && den_r > -den_i) || (den_r < den_i && den_r <= -den_i) {
            let r = den_i / den_r;
            let out_r = 1.0 / (den_r + r * den_i);
            let out_i = -r * out_r;
            (out_r, out_i)
        } else {
            let r = den_r / den_i;
            let out_i = -1.0 / (den_i + r * den_r);
            let out_r = -r * out_i;
            (out_r, out_i)
        }
    }

    fn element_mag(&self, idx: usize) -> f64 {
        let elem = &self.elements_arena[idx];
        if self.complex {
            elem.real.abs() + elem.imag.abs()
        } else {
            elem.real.abs()
        }
    }

    fn find_largest_in_col(&self, start_idx: Option<usize>) -> f64 {
        let mut largest = 0.0;
        let mut curr = start_idx;
        while let Some(idx) = curr {
            let mag = self.element_mag(idx);
            if mag > largest {
                largest = mag;
            }
            curr = self.elements_arena[idx].next_in_col;
        }
        largest
    }

    fn find_biggest_in_col_exclude(&self, p_element: usize, step: usize) -> f64 {
        let row = self.elements_arena[p_element].row;
        let col = self.elements_arena[p_element].col;

        let mut curr = self.first_in_col[col];
        while let Some(idx) = curr {
            if self.elements_arena[idx].row >= step {
                break;
            }
            curr = self.elements_arena[idx].next_in_col;
        }

        let mut largest = 0.0;
        if let Some(idx) = curr {
            if self.elements_arena[idx].row != row {
                largest = self.element_mag(idx);
            }
            curr = self.elements_arena[idx].next_in_col;
        }

        while let Some(idx) = curr {
            let mag = self.element_mag(idx);
            if mag > largest && self.elements_arena[idx].row != row {
                largest = mag;
            }
            curr = self.elements_arena[idx].next_in_col;
        }

        largest
    }

    fn real_row_col_elimination(&mut self, pivot_idx: usize) {
        if self.elements_arena[pivot_idx].real.abs() == 0.0 {
            self.matrix_is_singular(self.elements_arena[pivot_idx].row);
            return;
        }

        let pivot_val = 1.0 / self.elements_arena[pivot_idx].real;
        self.elements_arena[pivot_idx].real = pivot_val;

        let mut p_upper = self.elements_arena[pivot_idx].next_in_row;
        while let Some(u_idx) = p_upper {
            self.elements_arena[u_idx].real *= pivot_val;
            let u_real = self.elements_arena[u_idx].real;
            let u_col = self.elements_arena[u_idx].col;

            let mut p_sub = self.elements_arena[u_idx].next_in_col;
            let mut p_lower = self.elements_arena[pivot_idx].next_in_col;

            while let Some(l_idx) = p_lower {
                let row = self.elements_arena[l_idx].row;
                let l_real = self.elements_arena[l_idx].real;

                while let Some(sub_idx) = p_sub {
                    if self.elements_arena[sub_idx].row < row {
                        p_sub = self.elements_arena[sub_idx].next_in_col;
                    } else {
                        break;
                    }
                }

                let target_sub_idx = match p_sub {
                    Some(sub_idx) if self.elements_arena[sub_idx].row == row => sub_idx,
                    _ => match self.create_fillin(row, u_col) {
                        Some(fill_idx) => fill_idx,
                        None => {
                            self.error = 4; // spNO_MEMORY
                            return;
                        }
                    },
                };

                self.elements_arena[target_sub_idx].real -= u_real * l_real;
                p_sub = self.elements_arena[target_sub_idx].next_in_col;
                p_lower = self.elements_arena[l_idx].next_in_col;
            }
            p_upper = self.elements_arena[u_idx].next_in_row;
        }
    }

    fn complex_row_col_elimination(&mut self, pivot_idx: usize) {
        if self.element_mag(pivot_idx) == 0.0 {
            self.matrix_is_singular(self.elements_arena[pivot_idx].row);
            return;
        }

        let (piv_r, piv_i) = Self::cmplx_reciprocal(
            self.elements_arena[pivot_idx].real,
            self.elements_arena[pivot_idx].imag,
        );
        self.elements_arena[pivot_idx].real = piv_r;
        self.elements_arena[pivot_idx].imag = piv_i;

        let mut p_upper = self.elements_arena[pivot_idx].next_in_row;
        while let Some(u_idx) = p_upper {
            let ur0 = self.elements_arena[u_idx].real;
            let ui0 = self.elements_arena[u_idx].imag;
            let u_real = ur0 * piv_r - ui0 * piv_i;
            let u_imag = ur0 * piv_i + ui0 * piv_r;
            self.elements_arena[u_idx].real = u_real;
            self.elements_arena[u_idx].imag = u_imag;
            let u_col = self.elements_arena[u_idx].col;

            let mut p_sub = self.elements_arena[u_idx].next_in_col;
            let mut p_lower = self.elements_arena[pivot_idx].next_in_col;

            while let Some(l_idx) = p_lower {
                let row = self.elements_arena[l_idx].row;
                let l_real = self.elements_arena[l_idx].real;
                let l_imag = self.elements_arena[l_idx].imag;

                while let Some(sub_idx) = p_sub {
                    if self.elements_arena[sub_idx].row < row {
                        p_sub = self.elements_arena[sub_idx].next_in_col;
                    } else {
                        break;
                    }
                }

                let target_sub_idx = match p_sub {
                    Some(sub_idx) if self.elements_arena[sub_idx].row == row => sub_idx,
                    _ => match self.create_fillin(row, u_col) {
                        Some(fill_idx) => fill_idx,
                        None => {
                            self.error = 4; // spNO_MEMORY
                            return;
                        }
                    },
                };

                self.elements_arena[target_sub_idx].real -= u_real * l_real - u_imag * l_imag;
                self.elements_arena[target_sub_idx].imag -= u_real * l_imag + u_imag * l_real;
                p_sub = self.elements_arena[target_sub_idx].next_in_col;
                p_lower = self.elements_arena[l_idx].next_in_col;
            }
            p_upper = self.elements_arena[u_idx].next_in_row;
        }
    }

    fn count_markowitz(&mut self, rhs: Option<&[f64]>, step: usize) {
        // Generate MarkowitzRow Count for each row
        for i in step..=self.size {
            let mut count = -1_i32;
            let mut p_element = self.first_in_row[i];

            while let Some(idx) = p_element {
                if self.elements_arena[idx].col >= step {
                    count += 1;
                }
                p_element = self.elements_arena[idx].next_in_row;
            }

            if let Some(rhs_vec) = rhs {
                let ext_row = self.int_to_ext_row_map[i];
                if ext_row < rhs_vec.len() && rhs_vec[ext_row] != 0.0 {
                    count += 1;
                }
            }
            self.markowitz_row[i] = count.max(0) as usize;
        }

        // Generate MarkowitzCol count for each column
        for i in step..=self.size {
            let mut count = -1_i32;
            let mut p_element = self.first_in_col[i];

            while let Some(idx) = p_element {
                if self.elements_arena[idx].row >= step {
                    count += 1;
                }
                p_element = self.elements_arena[idx].next_in_col;
            }
            self.markowitz_col[i] = count.max(0) as usize;
        }
    }

    fn markowitz_products(&mut self, step: usize) {
        self.singletons = 0;

        for i in step..=self.size {
            let row_cnt = self.markowitz_row[i];
            let col_cnt = self.markowitz_col[i];

            let product = row_cnt.saturating_mul(col_cnt);
            self.markowitz_prod[i] = product;

            if product == 0 {
                self.singletons += 1;
            }
        }
    }

    fn search_for_singleton(&mut self, step: usize) -> Option<usize> {
        // Inspect `step` first (since choosing Diag[step] requires no row/col interchanges),
        // then `size` down to `step + 1` (1:1 with SearchForSingleton in spfactor.c).
        let mut singletons_left = self.singletons;
        if singletons_left == 0 {
            return None;
        }
        self.singletons -= 1;

        let iter = std::iter::once(step).chain((step + 1..=self.size).rev());
        for i in iter {
            if singletons_left == 0 {
                break;
            }
            if self.markowitz_prod[i] != 0 {
                continue;
            }
            singletons_left -= 1;

            if let Some(chosen_pivot) = self.diag[i] {
                let pivot_mag = self.element_mag(chosen_pivot);
                if pivot_mag > self.abs_threshold
                    && pivot_mag
                        > self.rel_threshold * self.find_biggest_in_col_exclude(chosen_pivot, step)
                {
                    return Some(chosen_pivot);
                }
            } else if self.markowitz_col[i] == 0 {
                let mut curr = self.first_in_col[i];
                while let Some(idx) = curr {
                    if self.elements_arena[idx].row >= step {
                        break;
                    }
                    curr = self.elements_arena[idx].next_in_col;
                }
                if let Some(chosen_pivot) = curr {
                    let pivot_mag = self.element_mag(chosen_pivot);
                    if pivot_mag > self.abs_threshold
                        && pivot_mag
                            > self.rel_threshold
                                * self.find_biggest_in_col_exclude(chosen_pivot, step)
                    {
                        return Some(chosen_pivot);
                    }
                }
                if self.markowitz_row[i] == 0 {
                    let mut rcurr = self.first_in_row[i];
                    while let Some(idx) = rcurr {
                        if self.elements_arena[idx].col >= step {
                            break;
                        }
                        rcurr = self.elements_arena[idx].next_in_row;
                    }
                    if let Some(chosen_pivot) = rcurr {
                        let pivot_mag = self.element_mag(chosen_pivot);
                        if pivot_mag > self.abs_threshold
                            && pivot_mag
                                > self.rel_threshold
                                    * self.find_biggest_in_col_exclude(chosen_pivot, step)
                        {
                            return Some(chosen_pivot);
                        }
                    }
                }
            } else {
                let mut rcurr = self.first_in_row[i];
                while let Some(idx) = rcurr {
                    if self.elements_arena[idx].col >= step {
                        break;
                    }
                    rcurr = self.elements_arena[idx].next_in_row;
                }
                if let Some(chosen_pivot) = rcurr {
                    let pivot_mag = self.element_mag(chosen_pivot);
                    if pivot_mag > self.abs_threshold
                        && pivot_mag
                            > self.rel_threshold
                                * self.find_biggest_in_col_exclude(chosen_pivot, step)
                    {
                        return Some(chosen_pivot);
                    }
                }
            }
        }

        // All singletons were unacceptable; restore count
        self.singletons += 1;
        None
    }

    fn quickly_search_diagonal(&mut self, step: usize) -> Option<usize> {
        const MAX_MARKOWITZ_TIES: usize = 100;
        const TIES_MULTIPLIER: usize = 5;

        let mut min_markowitz_product = usize::MAX;
        let mut tied_elements: Vec<usize> = Vec::with_capacity(MAX_MARKOWITZ_TIES + 1);

        let iter = std::iter::once(step).chain((step + 1..=self.size).rev());
        for i in iter {
            let prod = self.markowitz_prod[i];
            if prod > min_markowitz_product {
                continue;
            }

            let p_diag = match self.diag[i] {
                Some(idx) => idx,
                None => continue,
            };

            let magnitude = self.element_mag(p_diag);
            if magnitude <= self.abs_threshold {
                continue;
            }

            if prod == 1 {
                // Case where only one element exists in row and column other than diagonal
                let mut p_other_in_row = self.elements_arena[p_diag].next_in_row;
                let mut p_other_in_col = self.elements_arena[p_diag].next_in_col;
                if p_other_in_row.is_none() && p_other_in_col.is_none() {
                    let mut r = self.first_in_row[i];
                    while let Some(idx) = r {
                        let c = self.elements_arena[idx].col;
                        if c >= step && c != i {
                            p_other_in_row = Some(idx);
                            break;
                        }
                        r = self.elements_arena[idx].next_in_row;
                    }
                    let mut c = self.first_in_col[i];
                    while let Some(idx) = c {
                        let rw = self.elements_arena[idx].row;
                        if rw >= step && rw != i {
                            p_other_in_col = Some(idx);
                            break;
                        }
                        c = self.elements_arena[idx].next_in_col;
                    }
                }

                if let (Some(or_idx), Some(oc_idx)) = (p_other_in_row, p_other_in_col) {
                    if self.elements_arena[or_idx].col == self.elements_arena[oc_idx].row {
                        let largest_off_diag =
                            self.element_mag(or_idx).max(self.element_mag(oc_idx));
                        if magnitude >= largest_off_diag {
                            return Some(p_diag);
                        }
                    }
                }
            }

            if prod < min_markowitz_product {
                tied_elements.clear();
                tied_elements.push(p_diag);
                min_markowitz_product = prod;
            } else if tied_elements.len() <= MAX_MARKOWITZ_TIES {
                tied_elements.push(p_diag);
                if tied_elements.len() - 1 >= min_markowitz_product.saturating_mul(TIES_MULTIPLIER)
                {
                    break;
                }
            }
        }

        if tied_elements.is_empty() {
            return None;
        }

        let mut chosen_pivot = None;
        let mut max_ratio = 1.0 / self.rel_threshold;

        for &p_diag in &tied_elements {
            let magnitude = self.element_mag(p_diag);
            let largest_in_col = self.find_biggest_in_col_exclude(p_diag, step);
            let ratio = largest_in_col / magnitude;
            if ratio < max_ratio {
                chosen_pivot = Some(p_diag);
                max_ratio = ratio;
            }
        }

        chosen_pivot
    }

    fn search_diagonal(&mut self, step: usize) -> Option<usize> {
        const TIES_MULTIPLIER: usize = 5;

        let mut chosen_pivot = None;
        let mut min_markowitz_product = usize::MAX;
        let mut number_of_ties = 0usize;
        let mut ratio_of_accepted = 0.0;

        let iter = std::iter::once(step).chain((step + 1..=self.size).rev());
        for i in iter {
            let prod = self.markowitz_prod[i];
            if prod > min_markowitz_product {
                continue;
            }

            let p_diag = match self.diag[i] {
                Some(idx) => idx,
                None => continue,
            };

            let magnitude = self.element_mag(p_diag);
            if magnitude <= self.abs_threshold {
                continue;
            }

            let largest_in_col = self.find_biggest_in_col_exclude(p_diag, step);
            if magnitude <= self.rel_threshold * largest_in_col {
                continue;
            }

            if prod < min_markowitz_product {
                chosen_pivot = Some(p_diag);
                min_markowitz_product = prod;
                ratio_of_accepted = largest_in_col / magnitude;
                number_of_ties = 0;
            } else {
                number_of_ties += 1;
                let ratio = largest_in_col / magnitude;
                if ratio < ratio_of_accepted {
                    chosen_pivot = Some(p_diag);
                    ratio_of_accepted = ratio;
                }
                if number_of_ties >= min_markowitz_product.saturating_mul(TIES_MULTIPLIER) {
                    return chosen_pivot;
                }
            }
        }

        chosen_pivot
    }

    fn search_for_pivot(&mut self, step: usize, diag_pivoting: bool) -> Option<usize> {
        let mut chosen_pivot;

        if self.singletons > 0 {
            chosen_pivot = self.search_for_singleton(step);
            if chosen_pivot.is_some() {
                self.pivot_selection_method = 's';
                return chosen_pivot;
            }
        }

        if diag_pivoting {
            chosen_pivot = self.quickly_search_diagonal(step);
            if chosen_pivot.is_some() {
                self.pivot_selection_method = 'q';
                return chosen_pivot;
            }

            chosen_pivot = self.search_diagonal(step);
            if chosen_pivot.is_some() {
                self.pivot_selection_method = 'd';
                return chosen_pivot;
            }
        }

        chosen_pivot = self.search_entire_matrix(step);
        self.pivot_selection_method = 'e';
        chosen_pivot
    }
}
