use crate::domain_1d::Domain1D;
use crate::utils_csv::{read_csv, write_csv};
use crate::utils_error::FVChemError;
use crate::variable_1d::Variable1D;
use std::collections::HashMap;
use std::sync::Arc;

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
    pub is_write: bool,
    pub write_step: usize,
    pub write_file: String,

    // non-constant input type
    pub is_constant: bool,
    pub value_func: Arc<dyn Fn(usize, f64, &[f64]) -> f64 + Send + Sync>,
    pub value_var: Vec<usize>,
}

impl Default for Scalar1D {
    fn default() -> Self {
        Scalar1D {
            scl1d_id: 0,
            dom1d_id: 0,
            cell_value: HashMap::new(),
            cell_value_prev: HashMap::new(),
            face_value: HashMap::new(),
            face_value_prev: HashMap::new(),
            is_write: false,
            write_step: 0,
            write_file: String::new(),
            is_constant: true,
            value_func: Arc::new(|_, _, _| 0.0),
            value_var: Vec::new(),
        }
    }
}

impl Scalar1D {
    pub fn new(
        scl1d_id: usize,
        dom1d: &Domain1D,
        value: f64,
    ) -> Result<Scalar1D, FVChemError> {
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
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

        // set input to constant
        let is_constant = true;
        let value_func = Arc::new(|_: usize, _: f64, _: &[f64]| 0.0);
        let value_var: Vec<usize> = Vec::new();

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
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
        scl1d_id: usize,
        dom1d: &Domain1D,
        var1d_all: &Vec<Variable1D>,
        value_func: Arc<dyn Fn(usize, f64, &[f64]) -> f64 + Send + Sync>,
        value_var: Vec<usize>,
    ) -> Result<Scalar1D, FVChemError> {
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
            let scl_val = (value_func)(0, x, var_val.as_slice());
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
            let scl_val = (value_func)(0, x, var_val.as_slice());
            face_value.insert(fid, scl_val);
        }
        let face_value_prev: HashMap<i32, f64> = face_value.clone();

        // output file data
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

        // set input to non-constant
        let is_constant = false;

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
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
        scl1d_id: usize,
        dom1d: &Domain1D,
        value_file: String,
    ) -> Result<Scalar1D, FVChemError> {
        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // read cell data
        let (_, cell_i32, cell_f64) = read_csv(
            value_file.clone() + "_cell.csv",
            "Scalar1D".to_string(),
            vec!["cid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
        )?;
        let mut cell_value_raw: HashMap<i32, f64> = HashMap::new();
        for i in 0..cell_i32[0].len() {
            cell_value_raw.insert(cell_i32[0][i], cell_f64[1][i]);
        }
        let mut cell_value: HashMap<i32, f64> = HashMap::new();
        for &cid in dom1d.cell_id.iter() {
            cell_value.insert(cid, cell_value_raw[&cid]);
        }
        let cell_value_prev: HashMap<i32, f64> = cell_value.clone();

        // read face data
        let (_, face_i32, face_f64) = read_csv(
            value_file.clone() + "_face.csv",
            "Scalar1D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
        )?;
        let mut face_value_raw: HashMap<i32, f64> = HashMap::new();
        for i in 0..face_i32[0].len() {
            face_value_raw.insert(face_i32[0][i], face_f64[1][i]);
        }
        let mut face_value: HashMap<i32, f64> = HashMap::new();
        for &fid in dom1d.face_id.iter() {
            face_value.insert(fid, face_value_raw[&fid]);
        }
        let face_value_prev: HashMap<i32, f64> = face_value.clone();
        
        // output file data
        let is_write = false;
        let write_step = 0;
        let write_file = String::new();

        // set input to constant
        let is_constant = true;
        let value_func = Arc::new(|_: usize, _: f64, _: &[f64]| 0.0);
        let value_var: Vec<usize> = Vec::new();

        // return
        Ok(Scalar1D {
            scl1d_id,
            dom1d_id,
            cell_value,
            cell_value_prev,
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

    pub fn set_write_steady(scl: &mut Scalar1D, write_file: String) {
        scl.is_write = true;
        scl.write_file = write_file;
        scl.write_step = 0;
    }

    pub fn set_write_transient(scl: &mut Scalar1D, write_file: String, write_step: usize) {
        scl.is_write = true;
        scl.write_file = write_file;
        scl.write_step = write_step;
    }

    pub fn update_iter(dom: &Domain1D, scl: &mut Scalar1D, var_all: &Vec<Variable1D>, ts: usize) {
        // update only if non-constant
        if scl.is_constant {
            return;
        }

        // iterate over cells
        for &cid in dom.cell_id.iter() {
            // get position and variable values
            let x = dom.cell_x[&cid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &scl.value_var {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var_all[*v].cell_value[&cid]);
            }

            // update scalar value
            let scl_val = (&scl.value_func)(0, x, var_val.as_slice());
            scl.cell_value.insert(cid, scl_val);
        }

        // iterate over faces
        for &fid in dom.face_id.iter() {
            // get position and variable values
            let x = dom.face_x[&fid];
            let mut var_val: Vec<f64> = Vec::new();
            for v in &scl.value_var {
                // var_all[*v] -> Variable1D at index *v
                var_val.push(var_all[*v].face_value[&fid]);
            }

            // update scalar value
            let scl_val = (&scl.value_func)(ts, x, var_val.as_slice());
            scl.face_value.insert(fid, scl_val);
        }

    }

    pub fn update_prev(scl: &mut Scalar1D) {
        // copy current values to previous
        scl.cell_value_prev = scl.cell_value.clone();
        scl.face_value_prev = scl.face_value.clone();
    }

    pub fn write_steady(dom: &Domain1D, scl: &Scalar1D) -> Result<(), FVChemError> {
        // output only if specified
        if !scl.is_write {
            return Ok(());
        }

        // write cell data
        let mut cell_x_vec: Vec<f64> = Vec::new();
        let mut cell_value_vec: Vec<f64> = Vec::new();
        for &cid in dom.cell_id.iter() {
            cell_x_vec.push(dom.cell_x[&cid]);
            cell_value_vec.push(scl.cell_value[&cid]);
        }
        write_csv(
            scl.write_file.clone() + "_cell.csv",
            "Scalar1D".to_string(),
            vec!["cid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&dom.cell_id],
            vec![&cell_x_vec, &cell_value_vec],
        )?;

        // write face data
        let mut face_x_vec: Vec<f64> = Vec::new();
        let mut face_value_vec: Vec<f64> = Vec::new();
        for &fid in dom.face_id.iter() {
            face_x_vec.push(dom.face_x[&fid]);
            face_value_vec.push(scl.face_value[&fid]);
        }
        write_csv(
            scl.write_file.clone() + "_face.csv",
            "Scalar1D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&dom.face_id],
            vec![&face_x_vec, &face_value_vec],
        )?;

        // return
        Ok(())

    }

    pub fn write_transient(dom: &Domain1D, scl: &Scalar1D, ts: usize) -> Result<(), FVChemError> {
        // output only if specified
        if !scl.is_write || ts % scl.write_step != 0 {
            return Ok(());
        }

        // write cell data
        let mut cell_x_vec: Vec<f64> = Vec::new();
        let mut cell_value_vec: Vec<f64> = Vec::new();
        for &cid in dom.cell_id.iter() {
            cell_x_vec.push(dom.cell_x[&cid]);
            cell_value_vec.push(scl.cell_value[&cid]);
        }
        write_csv(
            scl.write_file.clone() + "_cell_" + &ts.to_string() + ".csv",
            "Scalar1D".to_string(),
            vec!["cid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&dom.cell_id],
            vec![&cell_x_vec, &cell_value_vec],
        )?;

        // write face data
        let mut face_x_vec: Vec<f64> = Vec::new();
        let mut face_value_vec: Vec<f64> = Vec::new();
        for &fid in dom.face_id.iter() {
            face_x_vec.push(dom.face_x[&fid]);
            face_value_vec.push(scl.face_value[&fid]);
        }
        write_csv(
            scl.write_file.clone() + "_face_" + &ts.to_string() + ".csv",
            "Scalar1D".to_string(),
            vec!["fid".to_string(), "x".to_string(), "u".to_string()],
            vec![true, false, false],
            vec![&dom.face_id],
            vec![&face_x_vec, &face_value_vec],
        )?;

        // return
        Ok(())

    }
}
