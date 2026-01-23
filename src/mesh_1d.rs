use crate::utils_csv::{read_csv, write_csv};
use crate::utils_error::FVChemError;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Mesh1D {
    // cell data
    pub num_cell: usize,
    pub cell_id: Vec<i32>,
    pub cell_x: HashMap<i32, f64>,
    pub cell_dx: HashMap<i32, f64>,

    // face data
    pub num_face: usize,
    pub face_id: Vec<i32>,
    pub face_x: HashMap<i32, f64>,

    // cell to cell data
    pub cell_cell_id: HashMap<i32, Vec<i32>>,
    pub cell_cell_dist: HashMap<i32, Vec<f64>>,

    // cell to face data
    pub cell_face_id: HashMap<i32, Vec<i32>>,
    pub cell_face_dist: HashMap<i32, Vec<f64>>,
    pub cell_face_norm: HashMap<i32, Vec<f64>>,

    // face to cell data
    pub face_cell_id: HashMap<i32, (i32, i32)>,
    pub face_cell_dist: HashMap<i32, (f64, f64)>,

    // region data
    pub num_reg: usize,
    pub reg_id: Vec<usize>,
    pub reg_cell_id: HashMap<usize, Vec<i32>>,

    // boundary data
    pub num_bnd: usize,
    pub bnd_id: Vec<usize>,
    pub bnd_cell_id: HashMap<usize, i32>,
    pub bnd_loc: HashMap<usize, usize>,
}

impl Mesh1D {
    pub fn new(x_min: f64, x_max: f64, num_cell: usize) -> Result<Mesh1D, FVChemError> {
        // error checking
        if x_max <= x_min {
            return Err(FVChemError::InvalidBoundsX {
                caller: "Mesh1D".to_string(),
                x_min: x_min,
                x_max: x_max,
            });
        }
        if num_cell < 1 {
            return Err(FVChemError::InvalidCellCount {caller: "Mesh1D".to_string()});
        }

        // cell data
        let mut cell_id: Vec<i32> = Vec::new();
        let mut cell_x: HashMap<i32, f64> = HashMap::new();
        let mut cell_dx: HashMap<i32, f64> = HashMap::new();

        // compute cell data
        let dx = (x_max - x_min) / (num_cell as f64);
        for i in 0..num_cell {
            let cid = i as i32; // cell id: 0 to num_cell-1
            let x = x_min + (i as f64 + 0.5) * dx;
            cell_id.push(cid);
            cell_x.insert(cid, x);
            cell_dx.insert(cid, dx);
        }

        // face data
        let mut face_id: Vec<i32> = Vec::new();
        let mut face_x: HashMap<i32, f64> = HashMap::new();

        // compute face data
        let num_face = num_cell + 1;
        for i in 0..num_face {
            let fid = -(i as i32 + 1); // face id: -1 to -num_face
            let x = x_min + (i as f64) * dx;
            face_id.push(fid);
            face_x.insert(fid, x);
        }

        // cell to cell data
        let mut cell_cell_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_cell_dist: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to cell data
        for i in 1..num_cell - 1 {
            // interior cells
            let cid = i as i32;
            cell_cell_id.insert(cid, vec![cid - 1, cid + 1]);
            cell_cell_dist.insert(cid, vec![dx, dx]);
        }
        cell_cell_id.insert(0, vec![-1, 1]); // no neighbor cell -> default to face
        cell_cell_dist.insert(0, vec![0.5 * dx, dx]);
        cell_cell_id.insert(
            num_cell as i32 - 1,
            vec![(num_cell as i32 - 2), -(num_face as i32)],
        );
        cell_cell_dist.insert(num_cell as i32 - 1, vec![dx, 0.5 * dx]);

        // cell to face data
        let mut cell_face_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_face_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        let mut cell_face_norm: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to face data
        for i in 0..num_cell {
            let cid = i as i32;
            cell_face_id.insert(cid, vec![-(i as i32 + 1), -(i as i32 + 2)]);
            cell_face_dist.insert(cid, vec![0.5 * dx, 0.5 * dx]);
            cell_face_norm.insert(cid, vec![-1.0, 1.0]);
        }

        // face to cell data
        let mut face_cell_id: HashMap<i32, (i32, i32)> = HashMap::new();
        let mut face_cell_dist: HashMap<i32, (f64, f64)> = HashMap::new();

        // compute face to cell data
        for i in 1..num_face - 1 {
            // interior faces
            let fid = -(i as i32 + 1);
            face_cell_id.insert(fid, (i as i32 - 1, i as i32));
            face_cell_dist.insert(fid, (0.5 * dx, 0.5 * dx));
        }
        face_cell_id.insert(-1, (-1, 0)); // no neighbor cell -> default to face
        face_cell_dist.insert(-1, (0.0, 0.5 * dx));
        face_cell_id.insert(
            -(num_face as i32),
            ((num_cell as i32 - 1), -(num_face as i32)),
        );
        face_cell_dist.insert(-(num_face as i32), (0.5 * dx, 0.0));

        // region data
        let num_reg = 1;
        let reg_id: Vec<usize> = vec![0];
        let reg_cell_id: HashMap<usize, Vec<i32>> = HashMap::from([(0, cell_id.clone())]);

        // boundary data
        let num_bnd = 2;
        let bnd_id: Vec<usize> = vec![0, 1];
        let bnd_cell_id: HashMap<usize, i32> = HashMap::from([(0, 0), (1, (num_cell as i32 - 1))]);
        let bnd_loc: HashMap<usize, usize> = HashMap::from([(0, 0), (1, 1)]);

        // return
        Ok(Mesh1D {
            num_cell,
            cell_id,
            cell_x,
            cell_dx,
            num_face,
            face_id,
            face_x,
            cell_cell_id,
            cell_cell_dist,
            cell_face_id,
            cell_face_dist,
            cell_face_norm,
            face_cell_id,
            face_cell_dist,
            num_reg,
            reg_id,
            reg_cell_id,
            num_bnd,
            bnd_id,
            bnd_cell_id,
            bnd_loc,
        })
    }

    pub fn new_from_file(mesh_file: String) -> Result<Mesh1D, FVChemError> {
        // read cell data
        let (num_cell, cell_i32, cell_f64) = read_csv(
            mesh_file.clone() + "_cell.csv",
            "Mesh1D".to_string(),
            vec!["cid".to_string(), "x".to_string(), "dx".to_string()],
            vec![true, false, false],
        )?;
        let mut cell_id: Vec<i32> = Vec::new();
        let mut cell_x: HashMap<i32, f64> = HashMap::new();
        let mut cell_dx: HashMap<i32, f64> = HashMap::new();
        for i in 0..num_cell {
            let cid = cell_i32[0][i];
            cell_id.push(cid);
            cell_x.insert(cid, cell_f64[0][i]);
            cell_dx.insert(cid, cell_f64[1][i]);
        }

        // read face data
        let (num_face, face_i32, face_f64) = read_csv(
            mesh_file.clone() + "_face.csv",
            "Mesh1D".to_string(),
            vec!["fid".to_string(), "x".to_string()],
            vec![true, false],
        )?;
        let mut face_id: Vec<i32> = Vec::new();
        let mut face_x: HashMap<i32, f64> = HashMap::new();
        for i in 0..num_face {
            let fid = face_i32[0][i];
            face_id.push(fid);
            face_x.insert(fid, face_f64[0][i]);
        }

        // read cell to cell data
        let (_, cell_cell_i32, cell_cell_f64) = read_csv(
            mesh_file.clone() + "_cell_cell.csv",
            "Mesh1D".to_string(),
            vec![
                "cid".to_string(),
                "nid_0".to_string(),
                "nid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
            ],
            vec![true, true, true, false, false],
        )?;
        let mut cell_cell_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_cell_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        for i in 0..num_cell {
            let cid = cell_cell_i32[0][i];
            cell_cell_id.insert(cid, vec![cell_cell_i32[1][i], cell_cell_i32[2][i]]);
            cell_cell_dist.insert(cid, vec![cell_cell_f64[0][i], cell_cell_f64[1][i]]);
        }

        // read cell to face data
        let (_, cell_face_i32, cell_face_f64) = read_csv(
            mesh_file.clone() + "_cell_face.csv",
            "Mesh1D".to_string(),
            vec![
                "cid".to_string(),
                "fid_0".to_string(),
                "fid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
                "norm_0".to_string(),
                "norm_1".to_string(),
            ],
            vec![true, true, true, false, false, false, false],
        )?;
        let mut cell_face_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_face_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        let mut cell_face_norm: HashMap<i32, Vec<f64>> = HashMap::new();
        for i in 0..num_cell {
            let cid = cell_face_i32[0][i];
            cell_face_id.insert(cid, vec![cell_face_i32[1][i], cell_face_i32[2][i]]);
            cell_face_dist.insert(cid, vec![cell_face_f64[0][i], cell_face_f64[1][i]]);
            cell_face_norm.insert(cid, vec![cell_face_f64[2][i], cell_face_f64[3][i]]);
        }

        // read face to cell data
        let (_, face_cell_i32, face_cell_f64) = read_csv(
            mesh_file.clone() + "_face_cell.csv",
            "Mesh1D".to_string(),
            vec![
                "fid".to_string(),
                "cid_0".to_string(),
                "cid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
            ],
            vec![true, true, true, false, false],
        )?;
        let mut face_cell_id: HashMap<i32, (i32, i32)> = HashMap::new();
        let mut face_cell_dist: HashMap<i32, (f64, f64)> = HashMap::new();
        for i in 0..num_face {
            let fid = face_cell_i32[0][i];
            face_cell_id.insert(fid, (face_cell_i32[1][i], face_cell_i32[2][i]));
            face_cell_dist.insert(fid, (face_cell_f64[0][i], face_cell_f64[1][i]));
        }

        // region data
        let (_, reg_i32, _) = read_csv(
            mesh_file.clone() + "_reg.csv",
            "Mesh1D".to_string(),
            vec!["rid".to_string(), "cid".to_string()],
            vec![true, true],
        )?;
        let mut reg_set: HashSet<usize> = HashSet::new();
        let mut reg_cell_id: HashMap<usize, Vec<i32>> = HashMap::new();
        for i in 0..reg_i32[0].len() {
            let rid = reg_i32[0][i] as usize;
            let cid = reg_i32[1][i];
            reg_set.insert(rid);
            reg_cell_id.entry(rid).or_insert(Vec::new()).push(cid);
        }
        let reg_id: Vec<usize> = reg_set.into_iter().collect();
        let num_reg = reg_id.len();

        // boundary data
        let (_, bnd_i32, _) = read_csv(
            mesh_file.clone() + "_bnd.csv",
            "Mesh1D".to_string(),
            vec!["bid".to_string(), "cid".to_string(), "loc".to_string()],
            vec![true, true, true],
        )?;
        let mut bnd_set: HashSet<usize> = HashSet::new();
        let mut bnd_cell_id: HashMap<usize, i32> = HashMap::new();
        let mut bnd_loc: HashMap<usize, usize> = HashMap::new();
        for i in 0..bnd_i32[0].len() {
            let bid = bnd_i32[0][i] as usize;
            let cid = bnd_i32[1][i];
            let loc = bnd_i32[2][i] as usize;
            bnd_set.insert(bid);
            bnd_cell_id.insert(bid, cid);
            bnd_loc.insert(bid, loc);
        }
        let bnd_id: Vec<usize> = bnd_set.iter().cloned().collect();
        let num_bnd = bnd_id.len();

        // return
        Ok(Mesh1D {
            num_cell,
            cell_id,
            cell_x,
            cell_dx,
            num_face,
            face_id,
            face_x,
            cell_cell_id,
            cell_cell_dist,
            cell_face_id,
            cell_face_dist,
            cell_face_norm,
            face_cell_id,
            face_cell_dist,
            num_reg,
            reg_id,
            reg_cell_id,
            num_bnd,
            bnd_id,
            bnd_cell_id,
            bnd_loc,
        })
    }

    pub fn write(&self, write_file: String) -> Result<(), FVChemError> {
        // write cell data
        let mut cell_x_vec: Vec<f64> = Vec::new();
        let mut cell_dx_vec: Vec<f64> = Vec::new();
        for &cid in self.cell_id.iter() {
            cell_x_vec.push(self.cell_x[&cid]);
            cell_dx_vec.push(self.cell_dx[&cid]);
        }
        write_csv(
            write_file.clone() + "_cell.csv",
            "Mesh1D".to_string(),
            vec!["cid".to_string(), "x".to_string(), "dx".to_string()],
            vec![true, false, false],
            vec![&self.cell_id],
            vec![&cell_x_vec, &cell_dx_vec],
        )?;

        // write face data
        let mut face_x_vec: Vec<f64> = Vec::new();
        for &fid in self.face_id.iter() {
            face_x_vec.push(self.face_x[&fid]);
        }
        write_csv(
            write_file.clone() + "_face.csv",
            "Mesh1D".to_string(),
            vec!["fid".to_string(), "x".to_string()],
            vec![true, false],
            vec![&self.face_id],
            vec![&face_x_vec],
        )?;

        // write cell to cell data
        let mut nid_0_vec: Vec<i32> = Vec::new();
        let mut nid_1_vec: Vec<i32> = Vec::new();
        let mut dist_0_vec: Vec<f64> = Vec::new();
        let mut dist_1_vec: Vec<f64> = Vec::new();
        for &cid in self.cell_id.iter() {
            let nid_sub = &self.cell_cell_id[&cid];
            let dist_sub = &self.cell_cell_dist[&cid];
            nid_0_vec.push(nid_sub[0]);
            nid_1_vec.push(nid_sub[1]);
            dist_0_vec.push(dist_sub[0]);
            dist_1_vec.push(dist_sub[1]);
        }
        write_csv(
            write_file.clone() + "_cell_cell.csv",
            "Mesh1D".to_string(),
            vec![
                "cid".to_string(),
                "nid_0".to_string(),
                "nid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
            ],
            vec![true, true, true, false, false],
            vec![&self.cell_id, &nid_0_vec, &nid_1_vec],
            vec![&dist_0_vec, &dist_1_vec],
        )?;

        // write cell to face data
        let mut fid_0_vec: Vec<i32> = Vec::new();
        let mut fid_1_vec: Vec<i32> = Vec::new();
        let mut dist_0_vec: Vec<f64> = Vec::new();
        let mut dist_1_vec: Vec<f64> = Vec::new();
        let mut norm_0_vec: Vec<f64> = Vec::new();
        let mut norm_1_vec: Vec<f64> = Vec::new();
        for &cid in self.cell_id.iter() {
            let fid_sub = &self.cell_face_id[&cid];
            let dist_sub = &self.cell_face_dist[&cid];
            let norm_sub = &self.cell_face_norm[&cid];
            fid_0_vec.push(fid_sub[0]);
            fid_1_vec.push(fid_sub[1]);
            dist_0_vec.push(dist_sub[0]);
            dist_1_vec.push(dist_sub[1]);
            norm_0_vec.push(norm_sub[0]);
            norm_1_vec.push(norm_sub[1]);
        }
        write_csv(
            write_file.clone() + "_cell_face.csv",
            "Mesh1D".to_string(),
            vec![
                "cid".to_string(),
                "fid_0".to_string(),
                "fid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
                "norm_0".to_string(),
                "norm_1".to_string(),
            ],
            vec![true, true, true, false, false, false, false],
            vec![&self.cell_id, &fid_0_vec, &fid_1_vec],
            vec![&dist_0_vec, &dist_1_vec, &norm_0_vec, &norm_1_vec],
        )?;

        // write face to cell data
        let mut cid_0_vec: Vec<i32> = Vec::new();
        let mut cid_1_vec: Vec<i32> = Vec::new();
        let mut dist_0_vec: Vec<f64> = Vec::new();
        let mut dist_1_vec: Vec<f64> = Vec::new();
        for &fid in self.face_id.iter() {
            let cid_sub = &self.face_cell_id[&fid];
            let dist_sub = &self.face_cell_dist[&fid];
            cid_0_vec.push(cid_sub.0);
            cid_1_vec.push(cid_sub.1);
            dist_0_vec.push(dist_sub.0);
            dist_1_vec.push(dist_sub.1);
        }
        write_csv(
            write_file.clone() + "_face_cell.csv",
            "Mesh1D".to_string(),
            vec![
                "fid".to_string(),
                "cid_0".to_string(),
                "cid_1".to_string(),
                "dist_0".to_string(),
                "dist_1".to_string(),
            ],
            vec![true, true, true, false, false],
            vec![&self.face_id, &cid_0_vec, &cid_1_vec],
            vec![&dist_0_vec, &dist_1_vec],
        )?;

        // write region data
        let mut rid_vec: Vec<i32> = Vec::new();
        let mut cid_vec: Vec<i32> = Vec::new();
        for (&rid, cid_sub) in self.reg_cell_id.iter() {
            for &cid in cid_sub.iter() {
                rid_vec.push(rid as i32);
                cid_vec.push(cid);
            }
        }
        write_csv(
            write_file.clone() + "_reg.csv",
            "Mesh1D".to_string(),
            vec!["rid".to_string(), "cid".to_string()],
            vec![true, true],
            vec![&rid_vec, &cid_vec],
            vec![],
        )?;

        // write boundary data
        let mut bid_vec: Vec<i32> = Vec::new();
        let mut cid_vec: Vec<i32> = Vec::new();
        let mut loc_vec: Vec<i32> = Vec::new();
        for &bid in self.bnd_id.iter() {
            bid_vec.push(bid as i32);
            cid_vec.push(self.bnd_cell_id[&bid]);
            loc_vec.push(self.bnd_loc[&bid] as i32);
        }
        write_csv(
            write_file.clone() + "_bnd.csv",
            "Mesh1D".to_string(),
            vec![
                "bid".to_string(),
                "cid".to_string(),
                "loc".to_string(),
            ],
            vec![true, true, true],
            vec![&bid_vec, &cid_vec, &loc_vec],
            vec![],
        )?;

        // return
        Ok(())

    }
}
