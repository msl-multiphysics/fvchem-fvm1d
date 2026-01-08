use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main()
{

    // create problem
    let mut prob = Problem1D::new();

    // create mesh with varying cell sizes
    let len_mesh = 1.0;  // length of mesh
    let num_cell = 20;  // number of cells
    let cell_rat = 0.2;  // size ratio of last to first cell
    let mut x_face = vec![0.0];
    let q = (cell_rat as f64).powf(1.0/(num_cell as f64 - 1.0));
    let x0 = len_mesh * (q - 1.0) / (q.powi(num_cell as i32) - 1.0);
    for i in 0..num_cell { 
        let dx = x0 * q.powi(i as i32);
        let x_new = x_face[i] + dx;
        x_face.push(x_new);
    }
    let mesh = Mesh1D::new_from_face(x_face);

    // add domain
    let dom = Problem1D::add_dom1d(&mut prob, &mesh);
    let dom_l = Problem1D::add_dom0d(&mut prob, dom, 0, 0);
    let dom_r = Problem1D::add_dom0d(&mut prob, dom, 19, 1);

    // add properties
    create_dir_all("examples/output_mesh_fromface").unwrap();
    let c = Problem1D::add_var1d(&mut prob, dom, 1.0, "examples/output_mesh_fromface/c".to_string(), 0);
    let d = Problem1D::add_scl1d(&mut prob, dom, 0.1, "".to_string(), 0);
    let r = Problem1D::add_scl1d(&mut prob, dom, 2.0, "".to_string(), 0);
    let c_l = Problem1D::add_scl0d(&mut prob, dom_l, 1.0, "".to_string(), 0);
    let n_r = Problem1D::add_scl0d(&mut prob, dom_r, 0.5, "".to_string(), 0);

    // create steady diffusion solver
    let mut solver = SteadyDiff::new();
    solver.add_domain(dom, c, d, r);
    solver.add_boundary_concentration(dom_l, c_l);
    solver.add_boundary_flux(dom_r, n_r);
    solver.solve(&mut prob, 1000, 1e-6);
    
}
