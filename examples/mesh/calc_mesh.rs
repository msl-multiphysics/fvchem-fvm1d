use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // output directory
    create_dir_all("examples/output_mesh").unwrap();

    // create mesh
    let mesh = Mesh1D::new(0.0, 1.0, 20).unwrap();
    mesh.write("examples/output_mesh/mesh".to_string()).unwrap();
}
