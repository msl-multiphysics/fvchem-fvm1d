pub enum Input1D {
    Constant(f64),
    Function(fn(f64, f64, Vec<f64>) -> f64, Vec<usize>),
    File(String),
}
