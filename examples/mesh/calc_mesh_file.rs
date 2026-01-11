use fvchem_fvm1d::*;
use std::fs::create_dir_all;

fn main() {
    // output directory
    create_dir_all("examples/output_mesh_file").unwrap();
    
    // create mesh
    let mesh = Mesh1D::new_from_file("examples/input_mesh_file/mesh".to_string()).unwrap();
    mesh.write("examples/output_mesh_file/mesh".to_string()).unwrap();
}
