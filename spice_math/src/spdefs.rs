//! DATA STRUCTURE DEFINITIONS for Sparse Matrix operations.
//!
//! Author: Kenneth S. Kundert (Original C version), Migrated to Rust
//! Advising professor: Alberto Sangiovanni-Vincentelli, UC Berkeley
//!
//! This module contains common type definitions and data structures for the sparse
//! matrix routines. In the original C version, memory was manually managed via a
//! custom allocator (elements blocked via chunked allocation) and orthogonal linked
//! lists via raw pointers.
//!
//! # Safety and Architectural Decisions
//! For the Rust port, safety and correctness take precedence. Rather than using
//! raw pointers and custom bump allocators for matrix elements, we utilize an
//! Arena-based approach (a `Vec<MatrixElement>`). Orthogonal linked list pointers
//! are replaced by `Option<usize>` indices pointing into this Arena. This guarantees
//! 100% memory safety while preserving the algorithmic structure of the Berkeley code.

pub type RealNumber = f64;

/// COMPLEX NUMBER DATA STRUCTURE
///
/// Basic structure to hold a complex number for the sparse solver.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ComplexNumber {
    /// The real portion of the number.
    pub real: RealNumber,
    /// The imaginary portion of the number.
    pub imag: RealNumber,
}

/// MATRIX ELEMENT DATA STRUCTURE
///
/// Every nonzero element in the matrix is stored dynamically. In the C version,
/// this was a struct within a custom orthogonal linked list. In this Rust
/// implementation, the elements are stored in an Arena (a `Vec` inside `MatrixFrame`),
/// and the `next_in_row` / `next_in_col` fields are indices into that vector.
///
/// Elements have an implicit assumption of being part of a larger `MatrixFrame`.
#[derive(Debug, Clone, Default)]
pub struct MatrixElement {
    /// The real portion of the value of the element.
    pub real: RealNumber,
    /// The imaginary portion of the value of the element.
    pub imag: RealNumber,
    /// The row number of the element (internal 1-based or 0-based representation).
    pub row: usize,
    /// The column number of the element (internal 1-based or 0-based representation).
    pub col: usize,
    /// Index to the next element in the row to the right of this element.
    /// If this element is the last nonzero in the row, then it contains `None`.
    pub next_in_row: Option<usize>,
    /// Index to the next element in the column below this element.
    /// If this element is the last nonzero in the column, then it contains `None`.
    pub next_in_col: Option<usize>,
}

/// MATRIX FRAME DATA STRUCTURE
///
/// This structure contains all the vectors that support the orthogonal linked list
/// containing the matrix elements. Also included in this structure are numbers
/// and state that are used globally by the sparse matrix routines and are associated
/// with one particular matrix.
///
/// Memory management artifacts (like `AllocationListPtr`, `RecordsRemaining`,
/// `NextAvailElement`) from the C struct were omitted since Rust's `Vec` provides
/// dynamic memory growth safely.
#[derive(Debug, Clone)]
pub struct MatrixFrame {
    /// The absolute magnitude an element must have to be considered as a pivot candidate,
    /// except as a last resort.
    pub abs_threshold: RealNumber,

    /// The magnitude an element must have relative to others in its row
    /// to be considered as a pivot candidate, except as a last resort.
    pub rel_threshold: RealNumber,

    /// The flag which indicates whether the matrix is complex (true) or real.
    pub complex: bool,

    /// Number of rows and columns in the matrix. Does not change as matrix is factored.
    pub size: usize,

    /// The value of the largest external row or column number encountered.
    pub ext_size: usize,

    /// Indicates the number of internal rows and columns that have elements in them
    /// during translation.
    pub current_size: usize,

    /// The arena storing all matrix elements. Indices into this Vec represent
    /// pointers in the original C linked list.
    pub elements_arena: Vec<MatrixElement>,

    /// Array of indices (into `elements_arena`) pointing to the diagonal elements.
    pub diag: Vec<Option<usize>>,

    /// Array of indices pointing to the first nonzero element of the column corresponding to the index.
    pub first_in_col: Vec<Option<usize>>,

    /// Array of indices pointing to the first nonzero element of the row corresponding to the index.
    pub first_in_row: Vec<Option<usize>>,

    /// Array of flags, one for each column in matrix. If a flag is true then corresponding
    /// column in a complex matrix should be eliminated in `spFactor()` using direct addressing.
    pub do_cmplx_direct: Vec<bool>,

    /// Same as above, but for real matrices.
    pub do_real_direct: Vec<bool>,

    /// The number of original elements (total elements minus fill-ins) present in matrix.
    pub elements: usize,

    /// The number of fill-ins created during the factorization of the matrix.
    pub fillins: usize,

    /// Error status.
    pub error: i32,

    /// An array used to convert external columns number to internal column numbers.
    pub ext_to_int_col_map: Vec<usize>,
    /// An array used to convert external row numbers to internal row numbers.
    pub ext_to_int_row_map: Vec<usize>,
    /// An array used to convert internal column numbers to external column numbers.
    pub int_to_ext_col_map: Vec<usize>,
    /// An array used to convert internal row numbers to external row numbers.
    pub int_to_ext_row_map: Vec<usize>,

    /// Indicates if matrix has been factored.
    pub factored: bool,

    /// A constant that provides the sparse data structure with a signature.
    pub id: u32,

    /// Temporary storage used in the spSolve routines during forward and backward substitution.
    pub intermediate: Vec<RealNumber>,

    /// Flag indicating whether the Markowitz vectors and Intermediate vector have been created.
    pub internal_vectors_allocated: bool,

    /// Count of the non-zero elements excluding pivots for each row. Used for MarkowitzProd.
    pub markowitz_row: Vec<usize>,
    /// Count of the non-zero elements excluding pivots for each col. Used for MarkowitzProd.
    pub markowitz_col: Vec<usize>,
    /// The product of Markowitz row and col counts. Element with smallest product
    /// is the best pivot to use to maintain sparsity.
    pub markowitz_prod: Vec<usize>,

    /// The maximum number of off-diagonal elements in the rows of L (lower triangular matrix).
    /// Used when computing an estimate of the roundoff error.
    pub max_row_count_in_lower_tri: usize,

    /// Signifies that the matrix needs to be ordered or reordered.
    pub needs_ordering: bool,

    /// Flag that indicates the sum of row and column interchange counts is an odd number.
    /// Used when determining the sign of the determinant.
    pub number_of_interchanges_is_odd: bool,

    /// Flag indicates that columns have been partitioned into two groups (direct vs indirect).
    pub partitioned: bool,

    /// Column pivot was chosen from.
    pub pivots_original_col: usize,
    /// Row pivot was chosen from.
    pub pivots_original_row: usize,
    /// Selection method used.
    pub pivot_selection_method: char,

    /// Flag to track if we need to zero out imaginary terms when clearing the matrix.
    pub previous_matrix_was_complex: bool,

    /// Signifies that the matrix has been reordered.
    pub reordered: bool,

    /// A flag that indicates whether the row pointers exist.
    pub rows_linked: bool,

    /// Normally zero, but if matrix is singular, assigned the external column number of zero pivot.
    pub singular_col: usize,
    /// Normally zero, but if matrix is singular, assigned the external row number of zero pivot.
    pub singular_row: usize,

    /// Number of singletons available for pivoting.
    pub singletons: usize,

    /// A dummy MatrixElement used to store elements added to the 'ground' (row 0 or col 0).
    pub trash_can: MatrixElement,
}
