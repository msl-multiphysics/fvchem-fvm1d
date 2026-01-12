use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // create mesh
    create_dir_all("examples/output_steady_itrcont").unwrap();
    let mesh = Mesh1D::new_from_file("examples/input_mesh_2dom/mesh".to_string()).unwrap();

    // add domain
    let mut prob = Problem1D::new();
    let dom_a = prob.add_dom1d(&mesh, 0).unwrap();
    let dom_b = prob.add_dom1d(&mesh, 1).unwrap();
    let dom_al = prob.add_dom0d(&mesh, dom_a, 0).unwrap();
    let dom_ar = prob.add_dom0d(&mesh, dom_a, 1).unwrap();
    let dom_bl = prob.add_dom0d(&mesh, dom_b, 2).unwrap();
    let dom_br = prob.add_dom0d(&mesh, dom_b, 3).unwrap();

    // add properties
    let c_a = prob.add_var1d(dom_a, 1.0).unwrap();
    let d_a = prob.add_scl1d(dom_a, 0.2).unwrap();
    let r_a = prob.add_scl1d(dom_a, 2.0).unwrap();
    let c_b = prob.add_var1d(dom_b, 1.0).unwrap();
    let d_b = prob.add_scl1d(dom_b, 0.1).unwrap();
    let r_b = prob.add_scl1d(dom_b, 1.0).unwrap();
    let c_l = prob.add_scl0d(dom_al, 1.0).unwrap();
    let n_r = prob.add_scl0d(dom_br, 0.5).unwrap();
    prob.set_var1d_write_steady(c_a, "examples/output_steady_itrcont/ca".to_string());
    prob.set_var1d_write_steady(c_b, "examples/output_steady_itrcont/cb".to_string());

    // create solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom_a, c_a, d_a, r_a);
    solver.add_domain(dom_b, c_b, d_b, r_b);
    solver.add_boundary_concentration(dom_al, c_l);
    solver.add_boundary_flux(dom_br, n_r);
    solver.add_interface_continuity(dom_ar, dom_bl);
    solver.solve(&mut prob, 1000, 1e-6, 1.0).unwrap();
}
