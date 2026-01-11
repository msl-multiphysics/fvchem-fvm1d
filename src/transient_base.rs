use crate::error_1d::Error1D;
use crate::problem_1d::Problem1D;
use crate::scalar_0d::Scalar0D;
use crate::scalar_1d::Scalar1D;
use crate::variable_1d::Variable1D;
use faer::linalg::solvers::Solve;
use faer::prelude::Col;
use faer::sparse::{SparseColMat, Triplet};
use std::time::{Duration, Instant};

pub trait TransientBase {
    // matrix assembly - to be implemented in specific models
    fn assemble_matrix(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dt: f64,
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

    // transient solver
    fn solve(
        &self,
        prob: &mut Problem1D,
        dt: f64,
        num_ts: usize,
        max_iter: usize,
        tol_l2: f64,
        damping: f64,
    ) -> Result<(), Error1D> {
        let time_start = Instant::now();
        println!("Starting transient solver.");

        // error checking
        if dt <= 0.0 {
            return Err(Error1D::InvalidTimeStep {caller: "Problem1D".to_string(), dt})
        }
        if num_ts < 1 {
            return Err(Error1D::InvalidNumTimeSteps {caller: "Problem1D".to_string(), num_ts})
        }
        if max_iter < 2 {
            return Err(Error1D::InvalidMaxIter {caller: "Problem1D".to_string(), max_iter})
        }
        if damping <= 0.0 || damping > 1.0 {
            return Err(Error1D::InvalidDamping {caller: "Problem1D".to_string(), damping})
        }

        // initialize time measurement
        let mut time_assemble = Duration::ZERO;
        let mut time_solve = Duration::ZERO;
        let mut time_update = Duration::ZERO;
        let mut time_write = Duration::ZERO;

        // initialize solver vectors
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec: Col<f64> = Col::zeros(0);
        let mut x_iter_vec: Col<f64> = Col::zeros(0);
        let mut mat_size: usize = 0;

        // resize solver vectors
        self.resize_vector(prob, &mut b_vec, &mut x_iter_vec, &mut mat_size);

        // iterate over time steps
        for ts in 0..num_ts {
            let time_0 = Instant::now();

            // write scalars and variables
            self.write_scalar_variable(prob, ts);

            let time_1 = Instant::now();

            // iterate to convergence
            let mut iter = 0;
            while iter < max_iter {
                let time_i0 = Instant::now();

                // update scalars
                self.update_scalar_iter(prob, ts);

                let time_i1 = Instant::now();

                // clear matrices and vectors
                a_triplet.clear();
                b_vec.fill(0.0);
                self.assemble_matrix(prob, &mut a_triplet, &mut b_vec, dt);

                let time_i2 = Instant::now();

                // solve linear system
                let a_mat = SparseColMat::<usize, f64>::try_new_from_triplets(mat_size, mat_size, &a_triplet);
                let a_mat = match a_mat {
                    Ok(mat) => mat,
                    Err(_) => return Err(Error1D::FailedMatrixAssembly {caller: "TransientBase".to_string()}),
                };
                let lu = a_mat.sp_lu();
                let lu = match lu {
                    Ok(decomp) => decomp,
                    Err(_) => return Err(Error1D::FailedLUDecomposition {caller: "TransientBase".to_string()}),
                };
                let x_undamp_vec = lu.solve(&b_vec); // undamped solution vector
                let x_vec = (1.0 - damping) * &x_iter_vec + damping * x_undamp_vec; // damped solution vector

                let time_i3 = Instant::now();

                // update variables
                self.update_variable_iter(prob, &x_vec);

                let time_i4 = Instant::now();

                // check convergence
                let err_l2 = (&x_vec - &x_iter_vec).norm_l2();
                println!("Time Step: {ts}; Iteration: {iter}; Residual: {err_l2}.");
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
                return Err(Error1D::FailedConvergence {caller: "TransientBase".to_string(), max_iter})
            }

            let time_2 = Instant::now();

            // update for next time step
            self.update_prev(prob);

            let time_3 = Instant::now();

            // update time measurements
            time_write += time_1.duration_since(time_0);
            time_update += time_3.duration_since(time_2);
        }

        let time_end = Instant::now();
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

    fn update_scalar_iter(&self, prob: &mut Problem1D, ts: usize) {
        // iterate over scalars and update
        for scl0d in &mut prob.scl0d {
            // &prob.dom0d[scl0d.dom0d_id] -> Domain0D
            // &prob.var1d -> Vec<Variable1D>
            Scalar0D::update_iter(&prob.dom0d[scl0d.dom0d_id], scl0d, &prob.var1d, ts);
        }
        for scl1d in &mut prob.scl1d {
            Scalar1D::update_iter(&prob.dom1d[scl1d.dom1d_id], scl1d, &prob.var1d, ts);
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

    fn update_prev(&self, prob: &mut Problem1D) {
        // update previous values
        for scl0d in &mut prob.scl0d {
            Scalar0D::update_prev(scl0d);
        }
        for scl1d in &mut prob.scl1d {
            Scalar1D::update_prev(scl1d);
        }
        for var in &mut prob.var1d {
            Variable1D::update_prev(var);
        }
    }

    fn write_scalar_variable(&self, prob: &mut Problem1D, ts: usize) {
        // write scalars and variables
        for scl0d in &prob.scl0d {
            Scalar0D::write_transient(&prob.dom0d[scl0d.dom0d_id], scl0d, ts);
        }
        for scl1d in &prob.scl1d {
            Scalar1D::write_transient(&prob.dom1d[scl1d.dom1d_id], scl1d, ts);
        }
        for var in &prob.var1d {
            Variable1D::write_transient(&prob.dom1d[var.dom1d_id], var, ts);
        }
    }

}
