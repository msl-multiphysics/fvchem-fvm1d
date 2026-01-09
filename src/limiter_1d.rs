use crate::problem_1d::Problem1D;

pub enum LimiterType {
    Linear,
    Upwind,
    BarthJespersen,
    Venkatakrishnan,
}

pub fn calc_limiter(
    prob: &Problem1D,
    dom_id: usize,
    var_id: usize,
    cid: i32,
    limiter_type: &LimiterType,
) -> f64 {
    match limiter_type {
        LimiterType::Linear => 1.0,
        LimiterType::Upwind => 0.0,
        LimiterType::BarthJespersen => calc_limiter_barth_jespersen(prob, dom_id, var_id, cid),
        LimiterType::Venkatakrishnan => calc_limiter_venkatakrishnan(prob, dom_id, var_id, cid),
    }
}

fn calc_limiter_barth_jespersen(prob: &Problem1D, dom_id: usize, var_id: usize, cid: i32) -> f64 {
    // preliminary values
    let mut limiter: f64 = 1.0;
    let mut var_min = prob.var1d[var_id].cell_value[&cid];
    let mut var_max = prob.var1d[var_id].cell_value[&cid];

    // loop over faces to find min and max in neighboring cells
    for loc in 0..2 {
        let nid = prob.dom1d[dom_id].cell_cell_id[&cid][loc]; // defaults to face id if boundary
        let var_n = if nid >= 0 {
            prob.var1d[var_id].cell_value[&nid]
        } else {
            prob.var1d[var_id].face_value[&nid]
        };
        var_min = var_min.min(var_n);
        var_max = var_max.max(var_n);
    }

    // iterate over faces to calculate limiter
    for loc in 0..2 {
        // get cell and neighbor concentrations
        let nid = prob.dom1d[dom_id].cell_cell_id[&cid][loc]; // defaults to face id if boundary
        let var_c = prob.var1d[var_id].cell_value[&cid];
        let var_n = if nid >= 0 {
            prob.var1d[var_id].cell_value[&nid]
        } else {
            prob.var1d[var_id].face_value[&nid]
        };

        // get distances
        let dist_cn = prob.dom1d[dom_id].cell_cell_dist[&cid][loc];
        let dist_cf = prob.dom1d[dom_id].cell_face_dist[&cid][loc];

        // get preliminary face value
        let var_f = var_c + (var_n - var_c) * (dist_cf / dist_cn);

        // compute limiter for this face
        let delta_var = var_f - var_c + 1e-12; // avoid zero division
        let r = if var_f > var_c {
            (var_max - var_c) / delta_var
        } else if var_f < var_c {
            (var_min - var_c) / delta_var
        } else {
            1.0
        };
        let limiter_sub = r.min(1.0);

        // update overall limiter
        limiter = limiter.min(limiter_sub);
    }

    // return
    limiter
}

fn calc_limiter_venkatakrishnan(prob: &Problem1D, dom_id: usize, var_id: usize, cid: i32) -> f64 {
    // preliminary values
    let mut limiter: f64 = 1.0;
    let mut var_min = prob.var1d[var_id].cell_value[&cid];
    let mut var_max = prob.var1d[var_id].cell_value[&cid];

    // loop over faces to find min and max in neighboring cells
    for loc in 0..2 {
        let nid = prob.dom1d[dom_id].cell_cell_id[&cid][loc]; // defaults to face id if boundary
        let var_n = if nid >= 0 {
            prob.var1d[var_id].cell_value[&nid]
        } else {
            prob.var1d[var_id].face_value[&nid]
        };
        var_min = var_min.min(var_n);
        var_max = var_max.max(var_n);
    }

    // iterate over faces to calculate limiter
    for loc in 0..2 {
        // get cell and neighbor concentrations
        let nid = prob.dom1d[dom_id].cell_cell_id[&cid][loc]; // defaults to face id if boundary
        let var_c = prob.var1d[var_id].cell_value[&cid];
        let var_n = if nid >= 0 {
            prob.var1d[var_id].cell_value[&nid]
        } else {
            prob.var1d[var_id].face_value[&nid]
        };

        // get distances
        let dist_cn = prob.dom1d[dom_id].cell_cell_dist[&cid][loc];
        let dist_cf = prob.dom1d[dom_id].cell_face_dist[&cid][loc];

        // get preliminary face value
        let var_f = var_c + (var_n - var_c) * (dist_cf / dist_cn);

        // compute limiter for this face
        let delta_var = var_f - var_c + 1e-12; // avoid zero division
        let r = if var_f > var_c {
            (var_max - var_c) / delta_var
        } else if var_f < var_c {
            (var_min - var_c) / delta_var
        } else {
            1.0
        };
        let limiter_sub = (r * r + 2.0 * r) / (r * r + r + 2.0);

        // update overall limiter
        limiter = limiter.min(limiter_sub);
    }

    // return
    limiter
}
