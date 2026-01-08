use fvchem_fvm1d::*;
use std::fs::create_dir_all;

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
    create_dir_all("examples/output_steady_diff_bndmtrn").unwrap();
    let c = Problem1D::add_var1d(&mut prob, dom, 1.0, "examples/output_steady_diff_bndmtrn/c".to_string(), 0);
    let d = Problem1D::add_scl1d(&mut prob, dom, 0.1, "".to_string(), 0);
    let r = Problem1D::add_scl1d(&mut prob, dom, 2.0, "".to_string(), 0);
    let k_l = Problem1D::add_scl0d(&mut prob, dom_l, 0.2, "".to_string(), 0);
    let cext_l = Problem1D::add_scl0d(&mut prob, dom_l, 1.0, "".to_string(), 0);
    let c_r = Problem1D::add_scl0d(&mut prob, dom_r, 1.0, "".to_string(), 0);

    // create steady diffusion solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_masstransfer(dom_l, k_l, cext_l);
    solver.add_boundary_concentration(dom_r, c_r);
    solver.solve(&mut prob, 1000, 1e-6);
    
}
