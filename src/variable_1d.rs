use crate::domain_1d::Domain1D;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

pub struct Variable1D{
    // struct ids
    pub var1d_id: usize,
    pub dom1d_id: usize,

    // cell data
    pub cell_value : HashMap<i32, f64>,
    pub cell_value_prev : HashMap<i32, f64>,

    // face data
    pub face_value : HashMap<i32, f64>,
    pub face_value_prev : HashMap<i32, f64>,

    // cell/face id to solution vector index
    pub xid : HashMap<i32, usize>,

    // output file data
    pub is_output : bool,
    pub output_step : usize,
    pub output_file : String,
}

impl Variable1D{
    pub fn new(var1d_id: usize, dom1d: &Domain1D, value_init: f64, output_file: String, output_step: usize) -> Variable1D{
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // cell data
        let mut cell_value : HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter(){
            cell_value.insert(cid, value_init);
        }
        let cell_value_prev : HashMap<i32, f64> = cell_value.clone();

        // face data
        let mut face_value : HashMap<i32, f64> = HashMap::new();
        for &fid in dom1d.face_id.iter(){
            face_value.insert(fid, value_init);
        }
        let face_value_prev : HashMap<i32, f64> = face_value.clone();

        // placeholder for cell/face id to solution vector index
        let mut xid : HashMap<i32, usize> = HashMap::new();
        for &cid in dom1d.cell_id.iter(){
            xid.insert(cid, 0);
        }
        for &fid in dom1d.face_id.iter(){
            xid.insert(fid, 0);
        }

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // return
        Variable1D{
            var1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            xid,
            is_output,
            output_step,
            output_file,
        }
    }

    pub fn update_prev(var: &mut Variable1D) {

        // update previous values
        var.cell_value_prev = var.cell_value.clone();
        var.face_value_prev = var.face_value.clone();

    }

    pub fn write_steady(dom: &Domain1D, var: &Variable1D)
    {

        // output only if specified
        if !var.is_output {
            return;
        }

        // write cell data
        let mut file_cell = File::create(var.output_file.clone() + "_cell.csv").unwrap();
        writeln!(file_cell, "x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(file_cell, "{:.6},{:.6}", dom.cell_x[&cid], var.cell_value[&cid]).unwrap();
        }

        // write face data
        let mut file_face = File::create(var.output_file.clone() + "_face.csv").unwrap();
        writeln!(file_face, "x,u").unwrap();
        for &fid in dom.face_id.iter() {    
            writeln!(file_face, "{:.6},{:.6}", dom.face_x[&fid], var.face_value[&fid]).unwrap();
        }

    }

    pub fn write_transient(dom: &Domain1D, var: &Variable1D, ts: usize)
    {

        // output only if specified
        if !var.is_output || ts % var.output_step != 0 {
            return;
        }

        // write cell data
        let mut file_cell = File::create(var.output_file.clone() + "_cell_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_cell, "x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(file_cell, "{:.6},{:.6}", dom.cell_x[&cid], var.cell_value[&cid]).unwrap();
        }

        // write face data
        let mut file_face = File::create(var.output_file.clone() + "_face_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_face, "x,u").unwrap();
        for &fid in dom.face_id.iter() {    
            writeln!(file_face, "{:.6},{:.6}", dom.face_x[&fid], var.face_value[&fid]).unwrap();
        }

    }

}
