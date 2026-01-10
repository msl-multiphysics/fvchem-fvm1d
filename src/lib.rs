// physics modules
pub mod steady_diff;
pub mod steady_diffadv;
pub mod transient_diff;
pub mod transient_diffadv;

// physics re-exports
pub use crate::steady_diff::*;
pub use crate::steady_diffadv::*;
pub use crate::transient_diff::*;
pub use crate::transient_diffadv::*;

// main modules
pub mod domain_0d;
pub mod domain_1d;
pub mod error_1d;
pub mod limiter_1d;
pub mod mesh_1d;
pub mod problem_1d;
pub mod scalar_0d;
pub mod scalar_1d;
pub mod steady_base;
pub mod transient_base;
pub mod variable_1d;

// main re-exports
pub use crate::domain_0d::*;
pub use crate::domain_1d::*;
pub use crate::mesh_1d::*;
pub use crate::problem_1d::*;
pub use crate::scalar_0d::*;
pub use crate::scalar_1d::*;
pub use crate::steady_base::*;
pub use crate::transient_base::*;
pub use crate::variable_1d::*;
