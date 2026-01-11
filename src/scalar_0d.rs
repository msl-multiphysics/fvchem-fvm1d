use crate::domain_0d::Domain0D;
use crate::utils_csv::{read_csv, write_csv};
use crate::utils_error::Error1D;
use crate::variable_1d::Variable1D;

pub struct Scalar0D {
    // struct ids
    pub scl0d_id: usize,
    pub dom0d_id: usize,

    // face data
    pub face_value: f64,
    pub face_value_prev: f64,

    // output file data
    pub is_write: bool,
    pub write_step: usize,
    pub write_file: String,

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
    ) -> Result<Scalar0D, Error1D> {
        // get struct ids
        let dom0d_id = dom0d.dom0d_id;

        // face data
        let face_value = value;
        let face_value_prev = value;

        // output file data
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

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
            is_write,
            write_step,
            write_file,
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
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

        // set input to non-constant
        let is_constant = false;

        // return
        Ok(Scalar0D {
            scl0d_id,
            dom0d_id,
            face_value,
            face_value_prev,
            is_write,
            write_step,
            write_file,
            is_constant,
            value_func,
            value_var,
        })
    }

    pub fn new_from_file(
        scl0d_id: usize,
        dom0d: &Domain0D,
        value_file: String,
    ) -> Result<Scalar0D, Error1D> {
        // get struct ids
        let dom0d_id = dom0d.dom0d_id;

        // read face data
        let(_, _, face_f64) = read_csv(
            value_file.clone() + "_0d.csv",
            "Scalar0D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
        )?;
        let face_value = face_f64[1][0];
        let face_value_prev = face_value;

        // output file data
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

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
            is_write,
            write_step,
            write_file,
            is_constant,
            value_func,
            value_var,
        })
    }

    pub fn set_write_steady(scl: &mut Scalar0D, write_file: String) {
        scl.is_write = true;
        scl.write_file = write_file;
        scl.write_step = 0;
    }

    pub fn set_write_transient(scl: &mut Scalar0D, write_file: String, write_step: usize) {
        scl.is_write = true;
        scl.write_file = write_file;
        scl.write_step = write_step;
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

    pub fn write_steady(dom: &Domain0D, scl: &Scalar0D) -> Result<(), Error1D> {
        // output only if specified
        if !scl.is_write {
            return Ok(());
        }

        // write face data
        let face_id_vec = vec![dom.face_id];
        let face_x_vec = vec![dom.face_x];
        let face_value_vec = vec![scl.face_value];
        write_csv(
            scl.write_file.clone() + "_0d.csv",
            "Scalar0D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&face_id_vec],
            vec![&face_x_vec, &face_value_vec],
        )?;

        // return
        Ok(())
    }

    pub fn write_transient(dom: &Domain0D, scl: &Scalar0D, ts: usize) -> Result<(), Error1D> {
        // output only if specified
        if !scl.is_write || scl.write_step <= 0|| ts % scl.write_step != 0 {
            return Ok(());
        }

        // write face data
        let face_id_vec = vec![dom.face_id];
        let face_x_vec = vec![dom.face_x];
        let face_value_vec = vec![scl.face_value];
        write_csv(
            scl.write_file.clone() + "_0d_" + &ts.to_string() + ".csv",
            "Scalar0D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&face_id_vec],
            vec![&face_x_vec, &face_value_vec],
        )?;

        // return
        Ok(())
    }

}
