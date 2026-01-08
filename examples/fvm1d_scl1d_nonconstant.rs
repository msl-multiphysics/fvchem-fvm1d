use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn r_func(_t: f64, _x: f64, vars: Vec<f64>) -> f64 {
    let c = vars[0];
    let r = 0.5 * (c - 1.0);  // source term function
    r
}

fn main()
{

    // create problem
    let mut prob = Problem1D::new();

    // add domain
    let mesh = Mesh1D::new(0.0, 1.0, 20);
    let dom = Problem1D::add_dom1d(&mut prob, &mesh);
    let dom_l = Problem1D::add_dom0d(&mut prob, dom, 0, 0);
    let dom_r = Problem1D::add_dom0d(&mut prob, dom, 19, 1);

    // add properties
    create_dir_all("examples/output_scl1d_nonconstant").unwrap();
    let c = Problem1D::add_var1d(&mut prob, dom, 1.0, "examples/output_scl1d_nonconstant/c".to_string(), 0);
    let d = Problem1D::add_scl1d(&mut prob, dom, 0.1, "".to_string(), 0);
    let r = Problem1D::add_scl1d_nonconstant(&mut prob, dom, r_func, vec![c], "".to_string(), 0);
    let c_l = Problem1D::add_scl0d(&mut prob, dom_l, 1.0, "".to_string(), 0);
    let n_r = Problem1D::add_scl0d(&mut prob, dom_r, 0.5, "".to_string(), 0);

    // create steady diffusion solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_flux(dom_r, n_r);
    solver.solve(&mut prob, 10, 1e-6);
    
}
