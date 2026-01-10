use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create problem
    let mut prob = Problem1D::new();

    // add domain
    let mesh = Mesh1D::new(0.0, 1.0, 20).unwrap();
    let dom = Problem1D::add_dom1d(&mut prob, &mesh).unwrap();
    let dom_l = Problem1D::add_dom0d(&mut prob, dom, 0, 0).unwrap();
    let dom_r = Problem1D::add_dom0d(&mut prob, dom, 19, 1).unwrap();

    // add properties
    create_dir_all("examples/output_steady_diff_bndflux").unwrap();
    let c = Problem1D::add_var1d(
        &mut prob,
        dom,
        1.0,
        "examples/output_steady_diff_bndflux/c".to_string(),
        0,
    ).unwrap();
    let d = Problem1D::add_scl1d(&mut prob, dom, 0.1, "".to_string(), 0).unwrap();
    let r = Problem1D::add_scl1d(&mut prob, dom, 2.0, "".to_string(), 0).unwrap();
    let c_l = Problem1D::add_scl0d(&mut prob, dom_l, 1.0, "".to_string(), 0).unwrap();
    let n_r = Problem1D::add_scl0d(&mut prob, dom_r, 0.5, "".to_string(), 0).unwrap();
    
    // create steady diffusion solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_flux(dom_r, n_r);
    solver.solve(&mut prob, 1000, 1e-6, 1.0).unwrap();
}
