use crate::domain_1d::Domain1D;
use crate::error_1d::Error1D;
use crate::input_1d::Input1D;
use crate::variable_1d::Variable1D;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub struct Scalar1D {
    // struct ids
    pub scl1d_id: usize,
    pub dom1d_id: usize,

    // cell data
    pub cell_value: HashMap<i32, f64>,
    pub cell_value_prev: HashMap<i32, f64>,

    // face data
    pub face_value: HashMap<i32, f64>,
    pub face_value_prev: HashMap<i32, f64>,

    // output file data
    pub is_output: bool,
    pub output_step: usize,
    pub output_file: String,

    // input type
    pub input_type: Input1D,
}

impl Scalar1D {
    pub fn new(
        scl1d_id: usize,
        dom1d: &Domain1D,
        value: f64,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar1D, Error1D> {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

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

        // input type
        let input_type = Input1D::Constant(value);

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            input_type,
        })
    }

    pub fn new_from_function(
        scl1d_id: usize,
        dom1d: &Domain1D,
        var1d_all: &Vec<Variable1D>,
        value_func: fn(f64, f64, Vec<f64>) -> f64,
        value_var: Vec<usize>,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar1D, Error1D> {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // cell data
        let mut cell_value: HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter() {
            // get position and variable values
            let x = dom1d.cell_x[&cid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &value_var {
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
            for v in &value_var {
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

        // input type
        let input_type = Input1D::Function(value_func, value_var.clone());

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            input_type,
        })
    }

    pub fn new_from_file(
        scl1d_id: usize,
        dom1d: &Domain1D,
        value_file: String,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar1D, Error1D> {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // raw cell data
        let cell_file = File::open(value_file.clone() + "_cell.csv");
        let cell_file = match cell_file {
            Ok(f) => f,
            Err(_) => {
                return Err(Error1D::FileNotFound{caller: "Scalar1D".to_string(), file_path: value_file.clone() + "_cell.csv" });
            }
        };
        let cell_reader = BufReader::new(cell_file);
        let mut cell_value_raw: HashMap<i32, f64> = HashMap::new();
        for (i, line) in cell_reader.lines().enumerate() {
            let line = match line {
                Ok(l) => l,
                Err(_) => {
                    return Err(Error1D::FileReadError{caller: "Scalar1D".to_string(), file_path: value_file.clone() + "_cell.csv" });
                }
            };
            if i == 0 {
                continue;  // skip header
            }
            let parts: Vec<&str> = line.split(',').collect();
            let cid: i32 = parts[0].trim().parse().unwrap();
            let val: f64 = parts[2].trim().parse().unwrap();
            cell_value_raw.insert(cid, val);
        }

        // cell data
        let mut cell_value: HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter() {
            let val = match cell_value_raw.get(&cid) {
                Some(v) => *v,
                None => {return Err(Error1D::InvalidCellID{caller: "Scalar1D".to_string(), cid, parent: "file".to_string()});}
            };
            cell_value.insert(cid, val);
        }
        let cell_value_prev: HashMap<i32, f64> = cell_value.clone();

        // raw face data
        let face_file = File::open(value_file.clone() + "_face.csv");
        let face_file = match face_file {
            Ok(f) => f,
            Err(_) => {
                return Err(Error1D::FileNotFound{caller: "Scalar1D".to_string(), file_path: value_file.clone() + "_face.csv" });
            }
        };
        let face_reader = BufReader::new(face_file);
        let mut face_value_raw: HashMap<i32, f64> = HashMap::new();
        for (i, line) in face_reader.lines().enumerate() {
            let line = match line {
                Ok(l) => l,
                Err(_) => {
                    return Err(Error1D::FileReadError{caller: "Scalar1D".to_string(), file_path: value_file.clone() + "_face.csv" });
                }
            };
            if i == 0 {
                continue;  // skip header
            }
            let parts: Vec<&str> = line.split(',').collect();
            let fid: i32 = parts[0].trim().parse().unwrap();
            let val: f64 = parts[2].trim().parse().unwrap();
            face_value_raw.insert(fid, val);
        }

        // face data
        let mut face_value: HashMap<i32, f64> = HashMap::new();
        for &fid in dom1d.face_id.iter() {
            let val = match face_value_raw.get(&fid) {
                Some(v) => *v,
                None => {return Err(Error1D::InvalidFaceID{caller: "Scalar1D".to_string(), fid, parent: "file".to_string()});}
            };
            face_value.insert(fid, val);
        }
        let face_value_prev: HashMap<i32, f64> = face_value.clone();
        
        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // input type
        let input_type = Input1D::File(value_file.clone());

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            input_type,
        })
    }

    pub fn update_iter(dom: &Domain1D, scl: &mut Scalar1D, var_all: &Vec<Variable1D>) {
        // update based on input type
        match &scl.input_type {
            Input1D::Constant(_val) => {
                return;  // no update needed
            }
            Input1D::File(_file_name) => {
                return;  // no update needed
            }
            Input1D::Function(value_func, value_var) => {
                // iterate over cells
                for &cid in dom.cell_id.iter() {
                    // get position and variable values
                    let x = dom.cell_x[&cid];
                    let mut var_val: Vec<f64> = Vec::new();
                    for v in value_var {
                        // var_all[*v] -> Variable1D at index *v
                        var_val.push(var_all[*v].cell_value[&cid]);
                    }

                    // update scalar value
                    let scl_val = (value_func)(0.0, x, var_val.clone());
                    scl.cell_value.insert(cid, scl_val);
                }

                // iterate over faces
                for &fid in dom.face_id.iter() {
                    // get position and variable values
                    let x = dom.face_x[&fid];
                    let mut var_val: Vec<f64> = Vec::new();
                    for v in value_var {
                        // var_all[*v] -> Variable1D at index *v
                        var_val.push(var_all[*v].face_value[&fid]);
                    }

                    // update scalar value
                    let scl_val = (value_func)(0.0, x, var_val.clone());
                    scl.face_value.insert(fid, scl_val);
                }
            }
        }

    }

    pub fn update_prev(scl: &mut Scalar1D) {
        // copy current values to previous
        scl.cell_value_prev = scl.cell_value.clone();
        scl.face_value_prev = scl.face_value.clone();
    }

    pub fn write_steady(dom: &Domain1D, scl: &Scalar1D) {
        // output only if specified
        if !scl.is_output {
            return;
        }

        // write cell data
        let mut file_cell = File::create(scl.output_file.clone() + "_cell.csv").unwrap();
        writeln!(file_cell, "cid,x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(
                file_cell,
                "{},{:.6},{:.6}",
                cid, dom.cell_x[&cid], scl.cell_value[&cid]
            )
            .unwrap();
        }

        // write face data
        let mut file_face = File::create(scl.output_file.clone() + "_face.csv").unwrap();
        writeln!(file_face, "fid,x,u").unwrap();
        for &fid in dom.face_id.iter() {
            writeln!(
                file_face,
                "{},{:.6},{:.6}",
                fid, dom.face_x[&fid], scl.face_value[&fid]
            )
            .unwrap();
        }
    }

    pub fn write_transient(dom: &Domain1D, scl: &Scalar1D, ts: usize) {
        // output only if specified
        if !scl.is_output || ts % scl.output_step != 0 {
            return;
        }

        // write cell data
        let mut file_cell =
            File::create(scl.output_file.clone() + "_cell_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_cell, "x,u").unwrap();
        for &cid in dom.cell_id.iter() {
            writeln!(
                file_cell,
                "{},{:.6},{:.6}",
                cid, dom.cell_x[&cid], scl.cell_value[&cid]
            )
            .unwrap();
        }

        // write face data
        let mut file_face =
            File::create(scl.output_file.clone() + "_face_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_face, "x,u").unwrap();
        for &fid in dom.face_id.iter() {
            writeln!(
                file_face,
                "{},{:.6},{:.6}",
                fid, dom.face_x[&fid], scl.face_value[&fid]
            )
            .unwrap();
        }
    }
}
