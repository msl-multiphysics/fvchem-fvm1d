use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FVChemError {

    #[error("{caller}: Cell ID {cid} not found in {parent}.")]
    InvalidCellID {
        caller: String,
        cid: i32,
        parent: String,
    },

    #[error("{caller}: Face ID {fid} not found in {parent}.")]
    InvalidFaceID {
        caller: String,
        fid: i32,
        parent: String,
    },

    #[error("{caller}: Face ID {fid} is not on the boundary of {parent}.")]
    InvalidBoundaryFaceID {
        caller: String,
        fid: i32,
        parent: String,
    },

    #[error("{caller}: Local ID {loc} not found in cell.")]
    InvalidLocalID {
        caller: String,
        loc: usize,
    },

    #[error("{caller}: Minimum bound must be less than maximum. Got x_min = {x_min}; x_max = {x_max}.")]
    InvalidBoundsX {
        caller: String,
        x_min: f64,
        x_max: f64,
    },
    
    #[error("{caller}: Need at least one cell.")]
    InvalidCellCount {
        caller: String,
    },

    #[error("{caller}: Cell size must be greater than zero.")]
    InvalidCellSize {
        caller: String,
    },

    #[error("{caller}: Need at last two iterations. Got max_iter = {max_iter}.")]
    InvalidMaxIter {
        caller: String,
        max_iter: usize,
    },

    #[error("{caller}: Damping factor must be in range (0, 1]. Got damping = {damping}.")]
    InvalidDamping {
        caller: String,
        damping: f64,
    },

    #[error("{caller}: Failed to assemble A matrix.")]
    FailedMatrixAssembly {
        caller: String,
    },

    #[error("{caller}: Failed to perform LU decomposition.")]
    FailedLUDecomposition {
        caller: String,
    },
    
    #[error("{caller}: Failed to converge in {max_iter} iterations.")]
    FailedConvergence {
        caller: String,
        max_iter: usize,
    },

    #[error("{caller}: Time step size must be greater than zero. Got dt = {dt}.")]
    InvalidTimeStep {
        caller: String,
        dt: f64,
    },

    #[error("{caller}: Need at least one time step. Got num_ts = {num_ts}.")]
    InvalidNumTimeSteps {
        caller: String,
        num_ts: usize,
    },

    #[error("{caller}: File not found: {file_path}.")]
    FileNotFound {
        caller: String,
        file_path: String,
    },

    #[error("{caller}: Failed to read from file: {file_path}.")]
    FileReadError {
        caller: String,
        file_path: String,
    },

    #[error("{caller}: Failed to write to file: {file_path}.")]
    FileWriteError {
        caller: String,
        file_path: String,
    },

    #[error("{caller}: File header does not match expected pattern: {file_path}.")]
    InvalidFileHeader {
        caller: String,
        file_path: String,
    },

}
