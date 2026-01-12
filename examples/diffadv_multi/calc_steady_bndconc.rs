use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create mesh
    create_dir_all("examples/output_steady_bndconc").unwrap();
    let mesh = Mesh1D::new(0.0, 1.0, 50).unwrap();

    // add domain
    let mut prob = Problem1D::new();
    let dom = prob.add_dom1d(&mesh, 0).unwrap();
    let dom_l = prob.add_dom0d(&mesh, dom, 0).unwrap();
    let dom_r = prob.add_dom0d(&mesh, dom, 1).unwrap();

    // add properties
    let c1 = prob.add_var1d(dom, 1.0).unwrap();
    let c2 = prob.add_var1d(dom, 1.0).unwrap();
    let d11 = prob.add_scl1d(dom, 0.1).unwrap();
    let d22 = prob.add_scl1d(dom, 0.2).unwrap();
    let u1 = prob.add_scl1d(dom, 2.0).unwrap();
    let u2 = prob.add_scl1d(dom, 10.0).unwrap();
    let r1 = prob.add_scl1d(dom, 0.0).unwrap();
    let r2 = prob.add_scl1d(dom, 0.0).unwrap();
    let c1_l = prob.add_scl0d(dom_l, 1.0).unwrap();
    let c2_l = prob.add_scl0d(dom_l, 1.0).unwrap();
    let c1_r = prob.add_scl0d(dom_r, 0.5).unwrap();
    let c2_r = prob.add_scl0d(dom_r, 0.5).unwrap();
    prob.set_var1d_write_steady(c1, "examples/output_steady_bndconc/c1".to_string());
    prob.set_var1d_write_steady(c2, "examples/output_steady_bndconc/c2".to_string());

    // create steady diffusion solver
    let mut solver = SteadyDiffAdvMulti::new(2, LimiterType::Linear);
    solver.add_domain(dom, 0, c1, [(0, d11)].iter().cloned().collect(), u1, r1);
    solver.add_domain(dom, 1, c2, [(1, d22)].iter().cloned().collect(), u2, r2);
    solver.add_boundary_concentration(dom_l, 0, c1_l);
    solver.add_boundary_concentration(dom_l, 1, c2_l);
    solver.add_boundary_concentration(dom_r, 0, c1_r);
    solver.add_boundary_concentration(dom_r, 1, c2_r);
    solver.solve(&mut prob, 1000, 1e-6, 0.5).unwrap();
}
