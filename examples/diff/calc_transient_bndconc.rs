use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create mesh
    create_dir_all("examples/output_transient_bndconc").unwrap();
    let mesh = Mesh1D::new(0.0, 1.0, 20).unwrap();

    // add domain
    let mut prob = Problem1D::new();
    let dom = prob.add_dom1d(&mesh, 0).unwrap();
    let dom_l = prob.add_dom0d(&mesh, dom, 0).unwrap();
    let dom_r = prob.add_dom0d(&mesh, dom, 1).unwrap();

    // add properties
    let c = prob.add_var1d(dom, 1.0).unwrap();
    let d = prob.add_scl1d(dom, 0.1).unwrap();
    let r = prob.add_scl1d(dom, 2.0).unwrap();
    let c_l = prob.add_scl0d(dom_l, 1.0).unwrap();
    let c_r = prob.add_scl0d(dom_r, 0.5).unwrap();
    prob.set_var1d_write_transient(c, "examples/output_transient_bndconc/c".to_string(), 20);
    
    // create steady diffusion solver
    let mut solver = TransientDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_concentration(dom_r, c_r);
    solver.solve(&mut prob, 0.01, 101, 1000, 1e-6, 1.0).unwrap();
}
