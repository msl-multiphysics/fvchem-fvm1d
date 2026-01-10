use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create problem
    let mut prob = Problem1D::new();

    // create mesh
    let mesh = Mesh1D::new(0.0, 1.5, 30).unwrap();

    // create domains
    let mut cid_a: Vec<i32> = Vec::new();
    let mut cid_b: Vec<i32> = Vec::new();
    for i in 0..10 {
        cid_a.push(i);
    }
    for i in 10..30 {
        cid_b.push(i);
    }
    let dom_a = Problem1D::add_dom1d_from_subset(&mut prob, &mesh, cid_a).unwrap();
    let dom_b = Problem1D::add_dom1d_from_subset(&mut prob, &mesh, cid_b).unwrap();
    let dom_al = Problem1D::add_dom0d(&mut prob, dom_a, 0, 0).unwrap();
    let dom_ar = Problem1D::add_dom0d(&mut prob, dom_a, 9, 1).unwrap();
    let dom_bl = Problem1D::add_dom0d(&mut prob, dom_b, 10, 0).unwrap();
    let dom_br = Problem1D::add_dom0d(&mut prob, dom_b, 29, 1).unwrap();

    // add properties
    create_dir_all("examples/output_steady_diff_itrcont").unwrap();
    let c_a = Problem1D::add_var1d(
        &mut prob,
        dom_a,
        1.0,
        "examples/output_steady_diff_itrcont/ca".to_string(),
        0,
    ).unwrap();
    let d_a = Problem1D::add_scl1d(&mut prob, dom_a, 0.2, "".to_string(), 0).unwrap();
    let r_a = Problem1D::add_scl1d(&mut prob, dom_a, 2.0, "".to_string(), 0).unwrap();
    let c_b = Problem1D::add_var1d(
        &mut prob,
        dom_b,
        1.0,
        "examples/output_steady_diff_itrcont/cb".to_string(),
        0,
    ).unwrap();
    let d_b = Problem1D::add_scl1d(&mut prob, dom_b, 0.1, "".to_string(), 0).unwrap();
    let r_b = Problem1D::add_scl1d(&mut prob, dom_b, 1.0, "".to_string(), 0).unwrap();
    let c_l = Problem1D::add_scl0d(&mut prob, dom_al, 1.0, "".to_string(), 0).unwrap();
    let n_r = Problem1D::add_scl0d(&mut prob, dom_br, 0.5, "".to_string(), 0).unwrap();
    
    // create steady diffusion solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom_a, c_a, d_a, r_a);
    solver.add_domain(dom_b, c_b, d_b, r_b);
    solver.add_boundary_concentration(dom_al, c_l);
    solver.add_boundary_flux(dom_br, n_r);
    solver.add_interface_continuity(dom_ar, dom_bl);
    solver.solve(&mut prob, 1000, 1e-6, 1.0).unwrap();
}
