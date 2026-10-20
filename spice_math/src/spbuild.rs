use crate::spdefs::{MatrixElement, MatrixFrame};

/// Translation array expansion.
const EXPANDABLE: bool = true;
const TRANSLATE: bool = true;

impl MatrixFrame {
    /// CLEAR MATRIX (spClear)
    ///
    /// Sets every element of the matrix to zero and clears the error flag.
    pub fn clear(&mut self) {
        if self.previous_matrix_was_complex || self.complex {
            for i in (1..=self.size).rev() {
                let mut elem_idx = self.first_in_col[i];
                while let Some(idx) = elem_idx {
                    self.elements_arena[idx].real = 0.0;
                    self.elements_arena[idx].imag = 0.0;
                    elem_idx = self.elements_arena[idx].next_in_col;
                }
            }
        } else {
            for i in (1..=self.size).rev() {
                let mut elem_idx = self.first_in_col[i];
                while let Some(idx) = elem_idx {
                    self.elements_arena[idx].real = 0.0;
                    elem_idx = self.elements_arena[idx].next_in_col;
                }
            }
        }

        self.trash_can.real = 0.0;
        self.trash_can.imag = 0.0;

        self.error = 0;
        self.factored = false;
        self.singular_col = 0;
        self.singular_row = 0;
        self.previous_matrix_was_complex = self.complex;
    }

    /// SINGLE ELEMENT ADDITION TO MATRIX BY INDEX (spGetElement)
    ///
    /// Finds element [Row,Col] and returns its index in the `elements_arena`.
    /// If the element is not found, it is created and spliced into the matrix.
    /// Returns the index instead of a raw pointer.
    pub fn get_element(&mut self, row: usize, col: usize) -> Option<usize> {
        if row == 0 || col == 0 {
            // Because we can't return a pointer to `trash_can` easily while mutating `self`
            // and returning indices elsewhere, returning a specific marker index or
            // simply using the `trash_can` approach directly is needed. For now, we return `None`
            // to indicate it's the trash can, which the caller can handle.
            return None;
        }

        let (mut int_row, mut int_col) = (row, col);

        if TRANSLATE {
            self.translate(&mut int_row, &mut int_col);
        } else if EXPANDABLE && (int_row > self.size || int_col > self.size) {
            self.enlarge_matrix(std::cmp::max(int_row, int_col));
        }

        if int_row == int_col {
            if let Some(idx) = self.diag[int_row] {
                return Some(idx);
            }
        }

        // Element does not exist or does not reside along diagonal. Search column.
        self.find_element_in_col(int_row, int_col, true)
    }

    fn translate(&mut self, row: &mut usize, col: &mut usize) {
        let (ext_row, ext_col) = (*row, *col);

        if ext_row > self.ext_to_int_col_map.len() - 1
            || ext_col > self.ext_to_int_col_map.len() - 1
        {
            self.expand_translation_arrays(std::cmp::max(ext_row, ext_col));
        }

        if ext_row > self.ext_size || ext_col > self.ext_size {
            self.ext_size = std::cmp::max(ext_row, ext_col);
        }

        let mut int_row = self.ext_to_int_row_map[ext_row];
        if int_row == usize::MAX {
            self.current_size += 1;
            int_row = self.current_size;
            self.ext_to_int_row_map[ext_row] = int_row;
            self.ext_to_int_col_map[ext_row] = int_row;

            if int_row > self.size {
                self.enlarge_matrix(int_row);
            }

            self.int_to_ext_row_map[int_row] = ext_row;
            self.int_to_ext_col_map[int_row] = ext_row;
        }

        if int_row > self.size {
            self.size = int_row;
        }

        let mut int_col = self.ext_to_int_col_map[ext_col];
        if int_col == usize::MAX {
            self.current_size += 1;
            int_col = self.current_size;
            self.ext_to_int_row_map[ext_col] = int_col;
            self.ext_to_int_col_map[ext_col] = int_col;

            if int_col > self.size {
                self.enlarge_matrix(int_col);
            }

            self.int_to_ext_row_map[int_col] = ext_col;
            self.int_to_ext_col_map[int_col] = ext_col;
        }

        if int_col > self.size {
            self.size = int_col;
        }

        *row = int_row;
        *col = int_col;
    }

    fn expand_translation_arrays(&mut self, new_ext_size: usize) {
        let new_size = std::cmp::max(new_ext_size, self.ext_to_int_col_map.len() + 100);
        self.ext_to_int_col_map.resize(new_size + 1, usize::MAX);
        self.ext_to_int_row_map.resize(new_size + 1, usize::MAX);
    }

    fn enlarge_matrix(&mut self, new_size: usize) {
        let old_allocated_size = self.diag.len() - 1;
        let new_size = std::cmp::max(new_size, old_allocated_size + 100);

        self.diag.resize(new_size + 1, None);
        self.first_in_col.resize(new_size + 1, None);
        self.first_in_row.resize(new_size + 1, None);
        self.int_to_ext_col_map.resize(new_size + 1, 0);
        self.int_to_ext_row_map.resize(new_size + 1, 0);

        for i in old_allocated_size + 1..=new_size {
            self.int_to_ext_col_map[i] = i;
            self.int_to_ext_row_map[i] = i;
        }
    }

    /// FIND ELEMENT BY SEARCHING COLUMN (spcFindElementInCol)
    pub(crate) fn find_element_in_col(
        &mut self,
        row: usize,
        col: usize,
        create_if_missing: bool,
    ) -> Option<usize> {
        let mut curr_idx = self.first_in_col[col];
        let mut prev_idx = None;

        while let Some(idx) = curr_idx {
            if self.elements_arena[idx].row < row {
                prev_idx = Some(idx);
                curr_idx = self.elements_arena[idx].next_in_col;
            } else if self.elements_arena[idx].row == row {
                return Some(idx);
            } else {
                break;
            }
        }

        if create_if_missing {
            self.create_element(row, col, prev_idx, false)
        } else {
            None
        }
    }

    /// CREATE AND SPLICE ELEMENT INTO MATRIX (spcCreateElement)
    pub(crate) fn create_element(
        &mut self,
        row: usize,
        col: usize,
        last_in_col_idx: Option<usize>,
        fillin: bool,
    ) -> Option<usize> {
        let new_idx = self.elements_arena.len();

        let mut elem = MatrixElement {
            real: 0.0,
            imag: 0.0,
            row,
            col,
            next_in_row: None,
            next_in_col: None,
        };

        if fillin {
            self.fillins += 1;
        } else {
            self.elements += 1;
        }

        if self.rows_linked {
            if !fillin {
                self.needs_ordering = true;
            }

            if row == col {
                self.diag[row] = Some(new_idx);
            }

            // Splice into column
            if let Some(prev_idx) = last_in_col_idx {
                elem.next_in_col = self.elements_arena[prev_idx].next_in_col;
                self.elements_arena[prev_idx].next_in_col = Some(new_idx);
            } else {
                elem.next_in_col = self.first_in_col[col];
                self.first_in_col[col] = Some(new_idx);
            }

            // Splash into row
            let mut curr_row_idx = self.first_in_row[row];
            let mut prev_row_idx = None;
            while let Some(idx) = curr_row_idx {
                if self.elements_arena[idx].col < col {
                    prev_row_idx = Some(idx);
                    curr_row_idx = self.elements_arena[idx].next_in_row;
                } else {
                    break;
                }
            }

            if let Some(prev_idx) = prev_row_idx {
                elem.next_in_row = self.elements_arena[prev_idx].next_in_row;
                self.elements_arena[prev_idx].next_in_row = Some(new_idx);
            } else {
                elem.next_in_row = self.first_in_row[row];
                self.first_in_row[row] = Some(new_idx);
            }
        } else {
            // If rows are not linked yet, we just push without row-splicing.
            // Splicing will happen in spcLinkRows().
            if row == col {
                self.diag[row] = Some(new_idx);
            }
            if let Some(prev_idx) = last_in_col_idx {
                elem.next_in_col = self.elements_arena[prev_idx].next_in_col;
                self.elements_arena[prev_idx].next_in_col = Some(new_idx);
            } else {
                elem.next_in_col = self.first_in_col[col];
                self.first_in_col[col] = Some(new_idx);
            }
        }

        self.elements_arena.push(elem);
        Some(new_idx)
    }
}
