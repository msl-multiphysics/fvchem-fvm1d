use crate::domain_0d::Domain0D;
use crate::error_1d::Error1D;
use crate::variable_1d::Variable1D;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub struct Scalar0D {
    // struct ids
    pub scl0d_id: usize,
    pub dom0d_id: usize,

    // face data
    pub face_value: f64,
    pub face_value_prev: f64,

    // output file data
    pub is_output: bool,
    pub output_step: usize,
    pub output_file: String,

    // non-constant input type
    pub is_constant: bool,
    pub value_func: fn(usize, f64, Vec<f64>) -> f64,
    pub value_var: Vec<usize>,
}

impl Scalar0D {
    pub fn new(
        scl0d_id: usize,
        dom0d: &Domain0D,
        value: f64,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar0D, Error1D> {
        // get struct ids
        let dom0d_id = dom0d.dom0d_id;

        // face data
        let face_value = value;
        let face_value_prev = value;

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // set input to constant
        let is_constant = true;
        let value_func = |_: usize, _: f64, _: Vec<f64>| 0.0;
        let value_var: Vec<usize> = Vec::new();

        // return
        Ok(Scalar0D {
            scl0d_id,
            dom0d_id,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            is_constant,
            value_func,
            value_var,
        })
    }

    pub fn new_from_function(
        scl0d_id: usize,
        dom0d: &Domain0D,
        var1d_all: &Vec<Variable1D>,
        value_func: fn(usize, f64, Vec<f64>) -> f64,
        value_var: Vec<usize>,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar0D, Error1D> {
        // get struct ids
        let dom0d_id = dom0d.dom0d_id;

        // face data
        let fid = dom0d.face_id;
        let x = dom0d.face_x;
        let mut var_val: Vec<f64> = Vec::new();
        for v in &value_var {
            // var_all[*v] -> Variable1D at index *v
            var_val.push(var1d_all[*v].face_value[&fid]);
        }
        let face_value = (value_func)(0, x, var_val.clone());
        let face_value_prev = face_value;

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // set input to non-constant
        let is_constant = false;

        // return
        Ok(Scalar0D {
            scl0d_id,
            dom0d_id,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            is_constant,
            value_func,
            value_var,
        })
    }

    pub fn new_from_file(
        scl0d_id: usize,
        dom0d: &Domain0D,
        value_file: String,
        output_file: String,
        output_step: usize,
    ) -> Result<Scalar0D, Error1D> {
        // get struct ids
        let dom0d_id = dom0d.dom0d_id;

        // raw face data
        let face_file = File::open(value_file.clone() + "_0d.csv");
        let face_file = match face_file {
            Ok(f) => f,
            Err(_) => {
                return Err(Error1D::FileReadError {
                    caller: "Scalar0D::new_from_file".to_string(),
                    file_path: value_file.clone() + "_0d.csv",
                })
            }
        };
        let reader = BufReader::new(face_file);
        let mut lines = reader.lines();

        // skip header
        lines.next();

        // read face value
        let line = lines.next().unwrap().unwrap();
        let parts: Vec<&str> = line.split(',').collect();
        let face_value: f64 = parts[2].parse().unwrap();
        let face_value_prev = face_value;

        // output file data
        let is_output = !output_file.is_empty() || output_step > 0;

        // set input to constant
        let is_constant = true;
        let value_func = |_: usize, _: f64, _: Vec<f64>| 0.0;
        let value_var: Vec<usize> = Vec::new();

        // return
        Ok(Scalar0D {
            scl0d_id,
            dom0d_id,
            face_value,
            face_value_prev,
            is_output,
            output_step,
            output_file,
            is_constant,
            value_func,
            value_var,
        })
    }

    pub fn update_iter(dom: &Domain0D, scl: &mut Scalar0D, var_all: &Vec<Variable1D>, ts: usize) {
        // update only if non-constant
        if scl.is_constant {
            return;
        }

        // get position and variable values
        let fid = dom.face_id;
        let x = dom.face_x;
        let mut var_val: Vec<f64> = Vec::new();
        for v in &scl.value_var {
            // var_all[*v] -> Variable1D at index *v
            var_val.push(var_all[*v].face_value[&fid]);
        }

        // update scalar value
        let scl_val = (&scl.value_func)(ts, x, var_val.clone());
        scl.face_value = scl_val;
    }

    pub fn update_prev(scl: &mut Scalar0D) {
        // copy current value to previous
        scl.face_value_prev = scl.face_value;
    }

    pub fn write_steady(dom: &Domain0D, scl: &Scalar0D) {
        // output only if specified
        if !scl.is_output {
            return;
        }

        // write face data
        let mut file_face = File::create(scl.output_file.clone() + "_0d.csv").unwrap();
        writeln!(file_face, "fid,x,u").unwrap();
        writeln!(file_face, "{},{:.6},{:.6}", dom.face_id, dom.face_x, scl.face_value).unwrap();
    }

    pub fn write_transient(dom: &Domain0D, scl: &Scalar0D, ts: usize) {
        // output only if specified
        if !scl.is_output || scl.output_step <= 0|| ts % scl.output_step != 0 {
            return;
        }

        // write face data
        let mut file_face =
            File::create(scl.output_file.clone() + "_0d_" + &ts.to_string() + ".csv").unwrap();
        writeln!(file_face, "fid,x,u").unwrap();
        writeln!(file_face, "{},{:.6},{:.6}", dom.face_id, dom.face_x, scl.face_value).unwrap();
    }

}
