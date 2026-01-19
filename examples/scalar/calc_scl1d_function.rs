use fvchem_fvm1d::*;
use std::fs::create_dir_all;
use std::sync::Arc;

fn r_func(_t: usize, _x: f64, vars: &[f64]) -> f64 {
    let c = vars[0];
    let r = -1.0 * (c - 2.0); // source term function
    r
}

fn main() {
    // create mesh
    create_dir_all("examples/output_scl1d_function").unwrap();
    let mesh = Mesh1D::new(0.0, 1.0, 20).unwrap();

    // add domain
    let mut prob = Problem1D::new();
    let dom = prob.add_dom1d(&mesh, 0).unwrap();
    let dom_l = prob.add_dom0d(&mesh, dom, 0).unwrap();
    let dom_r = prob.add_dom0d(&mesh, dom, 1).unwrap();

    // add properties
    let c = prob.add_var1d(dom, 1.0).unwrap();
    let d = prob.add_scl1d(dom, 0.1).unwrap();
    let r = prob.add_scl1d_from_function(dom, Arc::new(r_func), vec![c]).unwrap();
    let c_l = prob.add_scl0d(dom_l, 1.0).unwrap();
    let n_r = prob.add_scl0d(dom_r, 0.5).unwrap();
    prob.set_var1d_write_steady(c, "examples/output_scl1d_function/c".to_string());
    prob.set_scl1d_write_steady(r, "examples/output_scl1d_function/r".to_string());

    // create solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_flux(dom_r, n_r);
    solver.solve(&mut prob, 1000, 1e-6, 0.1).unwrap();
}
