use crate::utils_error::FVChemError;
use crate::problem_1d::Problem1D;
use crate::scalar_0d::Scalar0D;
use crate::scalar_1d::Scalar1D;
use crate::variable_1d::Variable1D;
use faer::linalg::solvers::Solve;
use faer::prelude::Col;
use faer::sparse::{SparseColMat, Triplet};
use std::time::{Duration, Instant};

pub trait SteadyBase {
    // matrix assembly - to be implemented in specific models
    fn assemble_matrix(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
    );

    // add terms to A
    fn add_a(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        var_row: usize,
        var_col: usize,
        row: i32,
        col: i32,
        val: f64,
    ) {
        let xrow = prob.var1d[var_row].xid[&row];
        let xcol = prob.var1d[var_col].xid[&col];
        a_triplet.push(Triplet::new(xrow, xcol, val));
    }

    // add terms to b
    fn add_b(&self, prob: &Problem1D, b_vec: &mut Col<f64>, var_row: usize, row: i32, val: f64) {
        let xrow = prob.var1d[var_row].xid[&row];
        b_vec[xrow] += val;
    }

    // steady-state solver
    fn solve(&self, prob: &mut Problem1D, max_iter: usize, tol_l2: f64, damping: f64) -> Result<(), FVChemError> {
        let time_start = Instant::now();
        println!("Starting steady-state solver.");

        // error checking
        if max_iter < 2 {
            return Err(FVChemError::InvalidMaxIter {caller: "Problem1D".to_string(), max_iter})
        }
        if damping <= 0.0 || damping > 1.0 {
            return Err(FVChemError::InvalidDamping {caller: "Problem1D".to_string(), damping})
        }

        // initialize time measurement
        let mut time_assemble = Duration::ZERO;
        let mut time_solve = Duration::ZERO;
        let mut time_update = Duration::ZERO;

        // initialize solver vectors
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec: Col<f64> = Col::zeros(0);
        let mut x_iter_vec: Col<f64> = Col::zeros(0);
        let mut mat_size: usize = 0;

        // resize solver vectors
        self.resize_vector(prob, &mut b_vec, &mut x_iter_vec, &mut mat_size);

        // iterate to convergence
        let mut iter = 0;
        while iter < max_iter {
            let time_i0 = Instant::now();

            // update scalars
            self.update_scalar_iter(prob);

            let time_i1 = Instant::now();

            // clear matrices and vectors
            a_triplet.clear();
            b_vec.fill(0.0);
            self.assemble_matrix(prob, &mut a_triplet, &mut b_vec);

            let time_i2 = Instant::now();

            // solve linear system
            let a_mat = SparseColMat::<usize, f64>::try_new_from_triplets(mat_size, mat_size, &a_triplet);
            let a_mat = match a_mat {
                Ok(mat) => mat,
                Err(_) => return Err(FVChemError::FailedMatrixAssembly {caller: "SteadyBase".to_string()}),
            };
            let lu = a_mat.sp_lu();
            let lu = match lu {
                Ok(decomp) => decomp,
                Err(_) => return Err(FVChemError::FailedLUDecomposition {caller: "SteadyBase".to_string()}),
            };
            let x_undamp_vec = lu.solve(&b_vec); // undamped solution vector
            let x_vec = (1.0 - damping) * &x_iter_vec + damping * x_undamp_vec; // damped solution vector

            let time_i3 = Instant::now();

            // update variables
            self.update_variable_iter(prob, &x_vec);

            let time_i4 = Instant::now();

            // check convergence
            let err_l2 = (&x_vec - &x_iter_vec).norm_l2();
            println!("Iteration: {iter}; Residual: {err_l2}.");
            if err_l2 < tol_l2 {
                break;
            }

            // update iteration vector
            x_iter_vec = x_vec;

            // update time measurements
            time_assemble += time_i2.duration_since(time_i1);
            time_solve += time_i3.duration_since(time_i2);
            time_update += time_i1.duration_since(time_i0) + time_i4.duration_since(time_i3);
        
            // increment iteration
            iter += 1;
        }

        // check convergence
        if iter >= max_iter {
            return Err(FVChemError::FailedConvergence {caller: "SteadyBase".to_string(), max_iter})
        }

        let time_0 = Instant::now();

        // write scalars and variables
        self.write_scalar_variable(prob)?;

        let time_end = Instant::now();
        let time_write = time_end.duration_since(time_0);
        let time_total = time_end.duration_since(time_start);

        // output time measurement (Duration -> seconds as f64)
        println!("Solution completed!");
        println!("Total time: {:.6} s", time_total.as_secs_f64());
        println!("  Assembly time: {:.6} s", time_assemble.as_secs_f64());
        println!("  Solve time: {:.6} s", time_solve.as_secs_f64());
        println!("  Update time: {:.6} s", time_update.as_secs_f64());
        println!("  Write time: {:.6} s", time_write.as_secs_f64());

        // return
        Ok(())
    }

    fn resize_vector(
        &self,
        prob: &mut Problem1D,
        b_vec: &mut Col<f64>,
        x_iter_vec: &mut Col<f64>,
        mat_size: &mut usize,
    ) {
        // determine size of matrix
        *mat_size = 0 as usize;
        for var in &prob.var1d {
            let dom_id = var.dom1d_id;
            *mat_size += prob.dom1d[dom_id].num_cell as usize;
            *mat_size += prob.dom1d[dom_id].num_face as usize;
        }

        // resize vectors
        *b_vec = Col::zeros(*mat_size);
        *x_iter_vec = Col::zeros(*mat_size);

        // map variable indices to matrix indices
        let mut xid = 0 as usize;
        for var in &mut prob.var1d {
            let dom_id = var.dom1d_id;
            for cid in prob.dom1d[dom_id].cell_id.iter() {
                var.xid.insert(*cid, xid);
                xid += 1;
            }
            for fid in prob.dom1d[dom_id].face_id.iter() {
                var.xid.insert(*fid, xid);
                xid += 1;
            }
        }

        // load initial variable values into iteration vector
        for var in &prob.var1d {
            let dom_id = var.dom1d_id;
            for cid in prob.dom1d[dom_id].cell_id.iter() {
                let xid = var.xid[&cid];
                let value = var.cell_value[&cid];
                x_iter_vec[xid] = value;
            }
            for fid in prob.dom1d[dom_id].face_id.iter() {
                let xid = var.xid[&fid];
                let value = var.face_value[&fid];
                x_iter_vec[xid] = value;
            }
        }
    }

    fn update_scalar_iter(&self, prob: &mut Problem1D) {
        // iterate over scalars and update
        for scl0d in &mut prob.scl0d {
            // &prob.dom0d[scl0d.dom0d_id] -> Domain0D
            // &prob.var1d -> Vec<Variable1D>
            Scalar0D::update_iter(&prob.dom0d[scl0d.dom0d_id], scl0d, &prob.var1d, 0);
        }
        for scl1d in &mut prob.scl1d {
            Scalar1D::update_iter(&prob.dom1d[scl1d.dom1d_id], scl1d, &prob.var1d, 0);
        }
    }

    fn update_variable_iter(&self, prob: &mut Problem1D, x_vec: &Col<f64>) {
        // iterate over variables and update
        for var in &mut prob.var1d {
            let dom_id = var.dom1d_id;
            for cid in prob.dom1d[dom_id].cell_id.iter() {
                let xid = var.xid[&cid];
                let value = x_vec[xid];
                var.cell_value.insert(*cid, value);
            }
            for fid in prob.dom1d[dom_id].face_id.iter() {
                let xid = var.xid[&fid];
                let value = x_vec[xid];
                var.face_value.insert(*fid, value);
            }
        }
    }

    fn write_scalar_variable(&self, prob: &mut Problem1D) -> Result<(), FVChemError> {
        // write scalars and variables
        for scl0d in &prob.scl0d {
            Scalar0D::write_steady(&prob.dom0d[scl0d.dom0d_id], scl0d)?;
        }
        for scl1d in &prob.scl1d {
            Scalar1D::write_steady(&prob.dom1d[scl1d.dom1d_id], scl1d)?;
        }
        for var in &prob.var1d {
            Variable1D::write_steady(&prob.dom1d[var.dom1d_id], var)?;
        }
        Ok(())
    }

}
