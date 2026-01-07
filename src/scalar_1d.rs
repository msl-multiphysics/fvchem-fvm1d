use crate::domain_1d::Domain1D;
use crate::variable_1d::Variable1D;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

pub struct Scalar1D {
    // struct ids
    pub scl1d_id: usize,
    pub dom1d_id: usize,
    pub var1d_id: Vec<usize>,

    // cell data
    pub cell_value: HashMap<i32, f64>,
    pub cell_value_prev : HashMap<i32, f64>,

    // face data
    pub face_value: HashMap<i32, f64>,
    pub face_value_prev : HashMap<i32, f64>,

    // output file data
    pub is_output: bool,
    pub output_step: usize,
    pub output_file: String,

    // non-constant data
    pub is_constant: bool,
    pub value_func: fn(f64, f64, Vec<f64>) -> f64,
}

impl Scalar1D {
    pub fn new(
        scl1d_id: usize,
        dom1d: &Domain1D,
        value: f64,
        output_file: String,
        output_step: usize,
    ) -> Scalar1D {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;
        let var1d_id: Vec<usize> = Vec::new();

        // cell data
        let mut cell_value: HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter() {
            cell_value.insert(cid, value);
        }
        let cell_value_prev: HashMap<i32, f64> = cell_value.clone();

        // face data
        let mut face_value: HashMap<i32, f64> = HashMap::new();
        for &fid in dom1d.face_id.iter() {
            face_value.insert(fid, value);
        }
        let face_value_prev: HashMap<i32, f64> = face_value.clone();

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // non-constant data
        let is_constant = true;
        let value_func = |_, _, _| 0.0;

        // return
        Scalar1D {
            scl1d_id,
            dom1d_id,
            var1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            is_constant,
            value_func,
        }
    }
    
    pub fn new_nonconstant(
        scl1d_id: usize,
        dom1d: &Domain1D,
        var1d_all: &Vec<Variable1D>,
        var1d_id: Vec<usize>,
        value_func: fn(f64, f64, Vec<f64>) -> f64,
        output_file: String,
        output_step: usize,
    ) -> Scalar1D {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // cell data
        let mut cell_value: HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter() {

            // get position and variable values
            let x = dom1d.cell_x[&cid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &var1d_id {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var1d_all[*v].cell_value[&cid]);
            }

            // update scalar value
            let scl_val = (value_func)(0.0, x, var_val.clone());
            cell_value.insert(cid, scl_val);
        }
        let cell_value_prev: HashMap<i32, f64> = cell_value.clone();

        // face data
        let mut face_value: HashMap<i32, f64> = HashMap::new();
        for &fid in dom1d.face_id.iter() {

            // get position and variable values
            let x = dom1d.face_x[&fid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &var1d_id {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var1d_all[*v].face_value[&fid]);
            }

            // update scalar value
            let scl_val = (value_func)(0.0, x, var_val.clone());
            face_value.insert(fid, scl_val);
        }
        let face_value_prev: HashMap<i32, f64> = face_value.clone();

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // non-constant data
        let is_constant = false;

        // return
        Scalar1D {
            scl1d_id,
            dom1d_id,
            var1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            is_constant,
            value_func,
        }
    }

    pub fn update_iter(dom: &Domain1D, scl: &mut Scalar1D, var_all: &Vec<Variable1D>) {
        
        // only update if non-constant
        if scl.is_constant {
            return;
        }

        // iterate over cells
        for &cid in dom.cell_id.iter() {

            // get position and variable values
            let x = dom.cell_x[&cid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &scl.var1d_id {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var_all[*v].cell_value[&cid]);
            }

            // update scalar value
            let scl_val = (scl.value_func)(0.0, x, var_val.clone());
            scl.cell_value.insert(cid, scl_val);

        }

        // iterate over faces
        for &fid in dom.face_id.iter() {

            // get position and variable values
            let x = dom.face_x[&fid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &scl.var1d_id {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var_all[*v].face_value[&fid]);
            }

            // update scalar value
            let scl_val = (scl.value_func)(0.0, x, var_val.clone());
            scl.face_value.insert(fid, scl_val);

        }

    }

    pub fn update_prev(scl: &mut Scalar1D) {

        // copy current values to previous
        scl.cell_value_prev = scl.cell_value.clone();
        scl.face_value_prev = scl.face_value.clone();

    }

    pub fn write_steady(dom: &Domain1D, scl: &Scalar1D)
    {

        // output only if specified
        if !scl.is_output {
            return;
        }

        // write cell data
        let mut file_cell = File::create(scl.output_file.clone() + "_cell.csv").unwrap();
        writeln!(file_cell, "x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(file_cell, "{:.6},{:.6}", dom.cell_x[&cid], scl.cell_value[&cid]).unwrap();
        }

        // write face data
        let mut file_face = File::create(scl.output_file.clone() + "_face.csv").unwrap();
        writeln!(file_face, "x,u").unwrap();
        for &fid in dom.face_id.iter() {    
            writeln!(file_face, "{:.6},{:.6}", dom.face_x[&fid], scl.face_value[&fid]).unwrap();
        }

    }

    pub fn write_transient(dom: &Domain1D, scl: &Scalar1D, ts: usize)
    {

        // output only if specified
        if !scl.is_output || ts % scl.output_step != 0 {
            return;
        }

        // write cell data
        let mut file_cell = File::create(scl.output_file.clone() + "_cell_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_cell, "x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(file_cell, "{:.6},{:.6}", dom.cell_x[&cid], scl.cell_value[&cid]).unwrap();
        }

        // write face data
        let mut file_face = File::create(scl.output_file.clone() + "_face_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_face, "x,u").unwrap();
        for &fid in dom.face_id.iter() {    
            writeln!(file_face, "{:.6},{:.6}", dom.face_x[&fid], scl.face_value[&fid]).unwrap();
        }

    }

}
