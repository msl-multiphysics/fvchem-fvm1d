use crate::utils_error::Error1D;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub fn write_csv(file_path: String, caller: String, header: Vec<String>, is_i32: Vec<bool>, data_i32: Vec<&Vec<i32>>, data_f64: Vec<&Vec<f64>>) -> Result<(), Error1D>
{
    // assume all columns have same number of rows
    // assume that is_i32 length matches data_i32 and data_f64 lengths
    // skip error checking for performance
    
    // open file
    let mut file = match File::create(&file_path) {
        Ok(f) => f,
        Err(_) => {
            return Err(Error1D::FileWriteError {caller, file_path});
        }
    };

    // write header
    let line_header = writeln!(file, "{}", header.join(","));
    if line_header.is_err() {
        return Err(Error1D::FileWriteError {caller, file_path});
    }

    // get number of rows and columns
    let num_col = is_i32.len();
    let num_row = if is_i32[0] {  // first column is i32
        data_i32[0].len()
    } else {  // first column is f64
        data_f64[0].len()
    };

    // iterate through rows
    for i in 0..num_row {

        // vector with string parts
        let mut parts: Vec<String> = Vec::with_capacity(num_col);

        // reset indices
        let mut idx_i32 = 0;
        let mut idx_f64 = 0;

        // iterate through columns
        for &is_i in is_i32.iter() {
            if is_i {  // i32 column
                let val = data_i32[idx_i32][i];
                parts.push(format!("{}", val));
                idx_i32 += 1;
            } else {  // f64 column
                let val = data_f64[idx_f64][i];
                parts.push(format!("{:.6}", val));
                idx_f64 += 1;
            }
        }

        // write line
        let line = writeln!(file, "{}", parts.join(","));
        if line.is_err() {
            return Err(Error1D::FileWriteError {caller: caller.clone(), file_path: file_path.clone()});
        }
    }

    // return
    Ok(())

}

pub fn read_csv(file_path: String, caller: String, header: Vec<String>, is_i32: Vec<bool>) -> Result<(usize, Vec<Vec<i32>>, Vec<Vec<f64>>), Error1D>
{
    // read file
    let file = File::open(&file_path);
    let file = match file {
        Ok(f) => f,
        Err(_) => {
            return Err(Error1D::FileReadError {caller, file_path});
        }
    };
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // check if header is consistent
    let line_header = match lines.next() {
        Some(Ok(l)) => l,
        _ => {
            return Err(Error1D::FileReadError {caller, file_path});
        }
    };
    if line_header.trim() != header.join(",").trim() {
        return Err(Error1D::InvalidFileHeader {caller, file_path});
    }

    // initialize data storage
    let mut num_i32 = 0;
    let mut num_f64 = 0;
    for &is_i in is_i32.iter() {
        if is_i {
            num_i32 += 1;
        } else {
            num_f64 += 1;
        }
    }
    let mut data_i32: Vec<Vec<i32>> = vec![Vec::new(); num_i32];
    let mut data_f64: Vec<Vec<f64>> = vec![Vec::new(); num_f64];

    // iterate through lines
    for line in lines {
        
        // parse line
        let line = match line {
            Ok(l) => l,
            Err(_) => {
                return Err(Error1D::FileReadError {caller: caller.clone(), file_path: file_path.clone()});
            }
        };
        let parts: Vec<&str> = line.split(',').collect();

        // reset indices
        let mut idx_i32 = 0;
        let mut idx_f64 = 0;

        // parse data
        for (j, &is_i) in is_i32.iter().enumerate() {
            if is_i {  // parse as i32
                let val: i32 = match parts[j].trim().parse() {
                    Ok(v) => v,
                    Err(_) => {
                        return Err(Error1D::FileReadError {caller: caller.clone(), file_path: file_path.clone()});
                    }
                };
                data_i32[idx_i32].push(val);
                idx_i32 += 1;
            } else {  // parse as f64
                let val: f64 = match parts[j].trim().parse() {
                    Ok(v) => v,
                    Err(_) => {
                        return Err(Error1D::FileReadError {caller: caller.clone(), file_path: file_path.clone()});
                    }
                };
                data_f64[idx_f64].push(val);
                idx_f64 += 1;
            }
        }
    }

    // get number of rows
    let num_row = if is_i32[0] {  // first column is i32
        data_i32[0].len()
    } else {  // first column is f64
        data_f64[0].len()
    };

    // return
    Ok((num_row, data_i32, data_f64))
}
