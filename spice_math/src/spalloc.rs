use crate::spdefs::{MatrixElement, MatrixFrame};

/// MACRO REPLACEMENTS
const EXPANDABLE: bool = true;
const TRANSLATE: bool = true;
const MINIMUM_ALLOCATED_SIZE: usize = 6;
const SPACE_FOR_ELEMENTS: usize = 6;
const SPACE_FOR_FILL_INS: usize = 4;
const DEFAULT_THRESHOLD: f64 = 1.0e-3;

impl Default for MatrixFrame {
    fn default() -> Self {
        Self::new(0, false)
    }
}

impl MatrixFrame {
    /// # MATRIX ALLOCATION
    ///
    /// Allocates and initializes the data structures associated with a matrix.
    /// This is equivalent to `spCreate` in the original C library. Because we
    /// use Rust's `Vec` for memory tracking, we do not need to mimic the complex
    /// `AllocationListPtr` logic.
    ///
    /// In C, this returned a pointer. In Rust, we return `Self`. Errors
    /// from `spNO_MEMORY` are naturally handled by Rust's global allocator panics,
    /// or we can just guarantee safe initialization.
    pub fn new(size: usize, complex: bool) -> Self {
        let allocated_size = std::cmp::max(size, MINIMUM_ALLOCATED_SIZE);

        let mut matrix = MatrixFrame {
            abs_threshold: 0.0,
            rel_threshold: DEFAULT_THRESHOLD,
            complex,
            size,
            ext_size: size,
            current_size: 0,

            // In the C version, these were linked lists of manually allocated blocks.
            // Here, we just use Vecs pre-allocated to the expected capacities.
            elements_arena: Vec::with_capacity(
                SPACE_FOR_ELEMENTS * allocated_size + SPACE_FOR_FILL_INS * allocated_size,
            ),

            diag: vec![None; allocated_size + 1],
            first_in_col: vec![None; allocated_size + 1],
            first_in_row: vec![None; allocated_size + 1],

            do_cmplx_direct: Vec::new(),
            do_real_direct: Vec::new(),

            elements: 0,
            fillins: 0,
            error: 0, // spOKAY

            ext_to_int_col_map: vec![0; allocated_size + 1],
            ext_to_int_row_map: vec![0; allocated_size + 1],
            int_to_ext_col_map: (0..=allocated_size).collect(),
            int_to_ext_row_map: (0..=allocated_size).collect(),

            factored: false,
            id: 0x772773,
            intermediate: Vec::new(),
            internal_vectors_allocated: false,

            markowitz_row: Vec::new(),
            markowitz_col: Vec::new(),
            markowitz_prod: Vec::new(),
            max_row_count_in_lower_tri: 0,
            needs_ordering: true,
            number_of_interchanges_is_odd: false,
            partitioned: false,
            pivots_original_col: 0,
            pivots_original_row: 0,
            pivot_selection_method: 'n',
            previous_matrix_was_complex: complex,
            reordered: false,
            rows_linked: false,
            singular_col: 0,
            singular_row: 0,
            singletons: 0,
            trash_can: MatrixElement::default(),
        };

        if TRANSLATE {
            for i in 1..=allocated_size {
                matrix.ext_to_int_col_map[i] = usize::MAX; // Represents -1
                matrix.ext_to_int_row_map[i] = usize::MAX; // Represents -1
            }
            matrix.ext_to_int_col_map[0] = 0;
            matrix.ext_to_int_row_map[0] = 0;
        }

        matrix
    }

    /// Instead of `spDestroy`, rust will automatically drop the `MatrixFrame`,
    /// cleanly deallocating the vectors without manual traversal.

    /// # MATRIX SIZE (spGetSize)
    ///
    /// Returns the size of the matrix. Either the internal or external size of
    /// the matrix is returned.
    pub fn get_size(&self, external: bool) -> usize {
        if external { self.ext_size } else { self.size }
    }

    /// WHERE IS MATRIX SINGULAR (spWhereSingular)
    ///
    /// Returns the row and column number where the matrix was detected as singular.
    pub fn where_singular(&self) -> (usize, usize) {
        // Assume spSINGULAR = 1, spZERO_DIAG = 2
        if self.error == 1 || self.error == 2 {
            (self.singular_row, self.singular_col)
        } else {
            (0, 0)
        }
    }
}
