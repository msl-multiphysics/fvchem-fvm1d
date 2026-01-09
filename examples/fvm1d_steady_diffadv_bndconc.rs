use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create problem
    let mut prob = Problem1D::new();

    // add domain
    let mesh = Mesh1D::new(0.0, 1.0, 50);
    let dom = Problem1D::add_dom1d(&mut prob, &mesh);
    let dom_l = Problem1D::add_dom0d(&mut prob, dom, 0, 0);
    let dom_r = Problem1D::add_dom0d(&mut prob, dom, 49, 1);

    // add properties
    create_dir_all("examples/output_steady_diffadv_bndconc").unwrap();
    let c = Problem1D::add_var1d(
        &mut prob,
        dom,
        1.0,
        "examples/output_steady_diffadv_bndconc/c".to_string(),
        0,
    );
    let d = Problem1D::add_scl1d(&mut prob, dom, 0.10, "".to_string(), 0);
    let u = Problem1D::add_scl1d(&mut prob, dom, 2.0, "".to_string(), 0);
    let r = Problem1D::add_scl1d(&mut prob, dom, 0.0, "".to_string(), 0);
    let c_l = Problem1D::add_scl0d(&mut prob, dom_l, 1.0, "".to_string(), 0);
    let c_r = Problem1D::add_scl0d(&mut prob, dom_r, 0.5, "".to_string(), 0);

    // create steady diffusion solver
    let mut solver = SteadyDiffAdv::new(limiter_1d::LimiterType::Linear);
    solver.add_domain(dom, c, d, u, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_concentration(dom_r, c_r);
    solver.solve(&mut prob, 1000, 1e-6, 1.0);
}
