use fvchem_fvm1d::*;
use std::{fs::create_dir_all, vec};

fn r1_func(_ts: usize, _x: f64, vars: Vec<f64>) -> f64 {
    // unpack variables
    let c1 = vars[0];
    let c2 = vars[1];

    // return reaction rate
    2.0 * c1 + 0.5 * c2
}

fn r2_func(_ts: usize, _x: f64, vars: Vec<f64>) -> f64 {
    // unpack variables
    let c1 = vars[0];
    let c2 = vars[1];

    // return reaction rate
    0.5 * c1 + 1.0 * c2
}

fn main() {
    // create mesh
    create_dir_all("examples/output_steady_itrcont").unwrap();
    let mesh = Mesh1D::new_from_file("examples/input_mesh_2dom/mesh".to_string()).unwrap();

    // domain A

    // add domain
    let mut prob = Problem1D::new();
    let dom_a = prob.add_dom1d(&mesh, 0).unwrap();
    let dom_al = prob.add_dom0d(&mesh, dom_a, 0).unwrap();
    let dom_ar = prob.add_dom0d(&mesh, dom_a, 1).unwrap();

    // add properties
    let c1_a = prob.add_var1d(dom_a, 1.0).unwrap();
    let c2_a = prob.add_var1d(dom_a, 1.0).unwrap();
    let d11_a = prob.add_scl1d(dom_a, 1.0).unwrap();
    let d12_a = prob.add_scl1d(dom_a, 0.5).unwrap();
    let d21_a = prob.add_scl1d(dom_a, 0.0).unwrap();
    let d22_a = prob.add_scl1d(dom_a, 2.0).unwrap();
    let r1_a = prob.add_scl1d_from_function(dom_a, r1_func, vec![c1_a, c2_a]).unwrap();
    let r2_a = prob.add_scl1d_from_function(dom_a, r2_func, vec![c1_a, c2_a]).unwrap();
    let c1_al = prob.add_scl0d(dom_al, 1.0).unwrap();
    let c2_al = prob.add_scl0d(dom_al, 0.5).unwrap();
    prob.set_var1d_write_steady(c1_a, "examples/output_steady_itrcont/c1_a".to_string());
    prob.set_var1d_write_steady(c2_a, "examples/output_steady_itrcont/c2_a".to_string());

    // domain B

    // add domain
    let dom_b = prob.add_dom1d(&mesh, 1).unwrap();
    let dom_bl = prob.add_dom0d(&mesh, dom_b, 2).unwrap();
    let dom_br = prob.add_dom0d(&mesh, dom_b, 3).unwrap();

    // add properties
    let c1_b = prob.add_var1d(dom_b, 1.0).unwrap();
    let c2_b = prob.add_var1d(dom_b, 1.0).unwrap();
    let d11_b = prob.add_scl1d(dom_b, 1.0).unwrap();
    let d12_b = prob.add_scl1d(dom_b, 0.5).unwrap();
    let d21_b = prob.add_scl1d(dom_b, 0.0).unwrap();
    let d22_b = prob.add_scl1d(dom_b, 2.0).unwrap();
    let r1_b = prob.add_scl1d_from_function(dom_b, r1_func, vec![c1_b, c2_b]).unwrap();
    let r2_b = prob.add_scl1d_from_function(dom_b, r2_func, vec![c1_b, c2_b]).unwrap();
    let c1_br = prob.add_scl0d(dom_br, 2.0).unwrap();
    let c2_br = prob.add_scl0d(dom_br, 1.0).unwrap();
    prob.set_var1d_write_steady(c1_b, "examples/output_steady_itrcont/c1_b".to_string());
    prob.set_var1d_write_steady(c2_b, "examples/output_steady_itrcont/c2_b".to_string());

    // create solver
    let mut solver = SteadyDiffMulti::new(2);
    solver.add_domain(dom_a, 0, c1_a, [(0, d11_a), (1, d12_a)].iter().cloned().collect(), r1_a);
    solver.add_domain(dom_a, 1, c2_a, [(0, d21_a), (1, d22_a)].iter().cloned().collect(), r2_a);
    solver.add_domain(dom_b, 0, c1_b, [(0, d11_b), (1, d12_b)].iter().cloned().collect(), r1_b);
    solver.add_domain(dom_b, 1, c2_b, [(0, d21_b), (1, d22_b)].iter().cloned().collect(), r2_b);
    solver.add_boundary_concentration(dom_al, 0, c1_al);
    solver.add_boundary_concentration(dom_al, 1, c2_al);
    solver.add_boundary_concentration(dom_br, 0, c1_br);
    solver.add_boundary_concentration(dom_br, 1, c2_br);
    solver.add_interface_continuity(dom_ar, dom_bl, 0);
    solver.add_interface_continuity(dom_ar, dom_bl, 1);
    solver.solve(&mut prob, 1000, 1e-6, 0.5).unwrap();

}
