use crate::error_1d::Error1D;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

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

    // group cell data
    pub num_group_cell: usize,
    pub group_cell: Vec<usize>,
    pub group_cell_id: HashMap<usize, Vec<i32>>,

    // group face data
    pub num_group_face: usize,
    pub group_face: Vec<usize>,
    pub group_face_id: HashMap<usize, Vec<i32>>,

}

impl Mesh1D {
    pub fn new(x_min: f64, x_max: f64, num_cell: usize) -> Result<Mesh1D, Error1D> {
        // error checking
        if x_max <= x_min {
            return Err(Error1D::InvalidBoundsX {
                caller: "Mesh1D".to_string(),
                x_min: x_min,
                x_max: x_max,
            });
        }
        if num_cell < 1 {
            return Err(Error1D::InvalidCellCount {caller: "Mesh1D".to_string()});
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

        // group cell data
        let num_group_cell = 1;
        let group_cell: Vec<usize> = vec![0];
        let group_cell_id: HashMap<usize, Vec<i32>> = HashMap::from([(0, cell_id.clone())]);

        // group face data
        let num_group_face = 2;
        let group_face: Vec<usize> = vec![0, 1];
        let group_face_id: HashMap<usize, Vec<i32>> = HashMap::from([
            (0, vec![-1]),
            (1, vec![-(num_face as i32)]),
        ]);

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
            num_group_cell,
            group_cell,
            group_cell_id,
            num_group_face,
            group_face,
            group_face_id,
        })
    }

    pub fn new_from_file(mesh_file: String) -> Result<Mesh1D, Error1D> {
        // read cell data
        let cell_file = File::open(mesh_file.clone() + "_cell.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_cell.csv",
            })?;
        let cell_reader = BufReader::new(cell_file);
        let mut cell_id: Vec<i32> = Vec::new();
        let mut cell_x: HashMap<i32, f64> = HashMap::new();
        let mut cell_dx: HashMap<i32, f64> = HashMap::new();
        for (i, line) in cell_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let cid: i32 = parts[0].parse().unwrap();
            let x: f64 = parts[1].parse().unwrap();
            let dx: f64 = parts[2].parse().unwrap();
            cell_id.push(cid);
            cell_x.insert(cid, x);
            cell_dx.insert(cid, dx);
        }
        let num_cell = cell_id.len();

        // read face data
        let face_file = File::open(mesh_file.clone() + "_face.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_face.csv",
            })?;
        let face_reader = BufReader::new(face_file);
        let mut face_id: Vec<i32> = Vec::new();
        let mut face_x: HashMap<i32, f64> = HashMap::new();
        for (i, line) in face_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let fid: i32 = parts[0].parse().unwrap();
            let x: f64 = parts[1].parse().unwrap();
            face_id.push(fid);
            face_x.insert(fid, x);
        }
        let num_face = face_id.len();

        // read cell to cell data
        let cell_cell_file = File::open(mesh_file.clone() + "_cell_cell.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_cell_cell.csv",
            })?;
        let cell_cell_reader = BufReader::new(cell_cell_file);
        let mut cell_cell_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_cell_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        for (i, line) in cell_cell_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let cid: i32 = parts[0].parse().unwrap();
            let nid_0: i32 = parts[1].parse().unwrap();
            let nid_1: i32 = parts[2].parse().unwrap();
            let dist_0: f64 = parts[3].parse().unwrap();
            let dist_1: f64 = parts[4].parse().unwrap();
            cell_cell_id.insert(cid, vec![nid_0, nid_1]);
            cell_cell_dist.insert(cid, vec![dist_0, dist_1]);
        }

        // read cell to face data
        let cell_face_file = File::open(mesh_file.clone() + "_cell_face.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_cell_face.csv",
            })?;
        let cell_face_reader = BufReader::new(cell_face_file);
        let mut cell_face_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_face_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        let mut cell_face_norm: HashMap<i32, Vec<f64>> = HashMap::new();
        for (i, line) in cell_face_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let cid: i32 = parts[0].parse().unwrap();
            let fid_0: i32 = parts[1].parse().unwrap();
            let fid_1: i32 = parts[2].parse().unwrap();
            let dist_0: f64 = parts[3].parse().unwrap();
            let dist_1: f64 = parts[4].parse().unwrap();
            let norm_0: f64 = parts[5].parse().unwrap();
            let norm_1: f64 = parts[6].parse().unwrap();
            cell_face_id.insert(cid, vec![fid_0, fid_1]);
            cell_face_dist.insert(cid, vec![dist_0, dist_1]);
            cell_face_norm.insert(cid, vec![norm_0, norm_1]);
        }

        // read face to cell data
        let face_cell_file = File::open(mesh_file.clone() + "_face_cell.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_face_cell.csv",
            })?;
        let face_cell_reader = BufReader::new(face_cell_file);
        let mut face_cell_id: HashMap<i32, (i32, i32)> = HashMap::new();
        let mut face_cell_dist: HashMap<i32, (f64, f64)> = HashMap::new();
        for (i, line) in face_cell_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let fid: i32 = parts[0].parse().unwrap();
            let cid_0: i32 = parts[1].parse().unwrap();
            let cid_1: i32 = parts[2].parse().unwrap();
            let dist_0: f64 = parts[3].parse().unwrap();
            let dist_1: f64 = parts[4].parse().unwrap();
            face_cell_id.insert(fid, (cid_0, cid_1));
            face_cell_dist.insert(fid, (dist_0, dist_1));
        }

        // read group cell data
        let group_cell_file = File::open(mesh_file.clone() + "_group_cell.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_group_cell.csv",
            })?;
        let group_cell_reader = BufReader::new(group_cell_file);
        let mut group_cell_set: HashSet<usize> = HashSet::new();
        let mut group_cell_id: HashMap<usize, Vec<i32>> = HashMap::new();
        for (i, line) in group_cell_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let gc: usize = parts[0].parse().unwrap();
            let cid: i32 = parts[1].parse().unwrap();
            group_cell_set.insert(gc);
            group_cell_id.entry(gc).or_insert(Vec::new()).push(cid);
        }
        let group_cell: Vec<usize> = group_cell_set.into_iter().collect();
        let num_group_cell = group_cell.len();

        // read group face data
        let group_face_file = File::open(mesh_file.clone() + "_group_face.csv")
            .map_err(|_| Error1D::FileReadError {
                caller: "Mesh1D".to_string(),
                file_path: mesh_file.clone() + "_group_face.csv",
            })?;
        let group_face_reader = BufReader::new(group_face_file);
        let mut group_face_set: HashSet<usize> = HashSet::new();
        let mut group_face_id: HashMap<usize, Vec<i32>> = HashMap::new();
        for (i, line) in group_face_reader.lines().enumerate() {
            if i == 0 {
                continue;  // skip header
            }
            let line = line.unwrap();
            let parts: Vec<&str> = line.trim().split(',').collect();
            let gf: usize = parts[0].parse().unwrap();
            let fid: i32 = parts[1].parse().unwrap();
            group_face_set.insert(gf);
            group_face_id.entry(gf).or_insert(Vec::new()).push(fid);
        }
        let group_face: Vec<usize> = group_face_set.into_iter().collect();
        let num_group_face = group_face.len();

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
            num_group_cell,
            group_cell,
            group_cell_id,
            num_group_face,
            group_face,
            group_face_id,
        })
    }

    pub fn write(&self, write_file: String) -> Result<(), Error1D> {
        // write cell data
        let path_cell = write_file.clone() + "_cell.csv";
        let mut file_cell = File::create(&path_cell).map_err(|_| Error1D::FileWriteError {
            caller: "Mesh1D".to_string(),
            file_path: path_cell.clone(),
        })?;
        writeln!(file_cell, "cid,x,dx").unwrap();
        for &cid in &self.cell_id {
            writeln!(file_cell, "{},{:.6},{:.6}", cid, self.cell_x[&cid], self.cell_dx[&cid]).unwrap();
        }
        
        // write face data
        let path_face = write_file.clone() + "_face.csv";
        let mut file_face = File::create(&path_face)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_face.clone(),
            })?;
        writeln!(file_face, "fid,x").unwrap();
        for &fid in &self.face_id {
            writeln!(file_face, "{},{:.6}", fid, self.face_x[&fid]).unwrap();
        }

        // write cell to cell data
        let path_cell_cell = write_file.clone() + "_cell_cell.csv";
        let mut file_cell_cell = File::create(&path_cell_cell)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_cell_cell.clone(),
            })?;
        writeln!(file_cell_cell, "cid,nid_0,nid_1,dist_0,dist_1").unwrap();
        for &cid in self.cell_id.iter() {
            let nid_vec = &self.cell_cell_id[&cid];
            let dist_vec = &self.cell_cell_dist[&cid];
            writeln!(file_cell_cell, "{},{},{},{:.6},{:.6}", cid, nid_vec[0], nid_vec[1], dist_vec[0], dist_vec[1]).unwrap();
        }

        // write cell to face data
        let path_cell_face = write_file.clone() + "_cell_face.csv";
        let mut file_cell_face = File::create(&path_cell_face)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_cell_face.clone(),
            })?;
        writeln!(file_cell_face, "cid,fid_0,fid_1,dist_0,dist_1,norm_0,norm_1").unwrap();
        for &cid in self.cell_id.iter() {
            let fid_vec = &self.cell_face_id[&cid];
            let dist_vec = &self.cell_face_dist[&cid];
            let norm_vec = &self.cell_face_norm[&cid];
            writeln!(file_cell_face, "{},{},{},{:.6},{:.6},{:.6},{:.6}", cid, fid_vec[0], fid_vec[1], dist_vec[0], dist_vec[1], norm_vec[0], norm_vec[1]).unwrap();
        }

        // write face to cell data
        let path_face_cell = write_file.clone() + "_face_cell.csv";
        let mut file_face_cell = File::create(&path_face_cell)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_face_cell.clone(),
            })?;
        writeln!(file_face_cell, "fid,cid_0,cid_1,dist_0,dist_1").unwrap();
        for &fid in self.face_id.iter() {
            let (cid_0, cid_1) = self.face_cell_id[&fid];
            let (dist_0, dist_1) = self.face_cell_dist[&fid];
            writeln!(file_face_cell, "{},{},{},{:.6},{:.6}", fid, cid_0, cid_1, dist_0, dist_1).unwrap();
        }

        // write group cell data
        let path_group_cell = write_file.clone() + "_group_cell.csv";
        let mut file_group_cell = File::create(&path_group_cell)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_group_cell.clone(),
            })?;
        writeln!(file_group_cell, "gc,cid").unwrap();
        for &gc in self.group_cell.iter() {
            let cid_vec = &self.group_cell_id[&gc];
            for &cid in cid_vec.iter() {
                writeln!(file_group_cell, "{},{}", gc, cid).unwrap();
            }
        };

        // write group face data
        let path_group_face = write_file.clone() + "_group_face.csv";
        let mut file_group_face = File::create(&path_group_face)
            .map_err(|_| Error1D::FileWriteError {
                caller: "Mesh1D".to_string(),
                file_path: path_group_face.clone(),
            })?;
        writeln!(file_group_face, "gf,fid").unwrap();
        for &gf in self.group_face.iter() {
            let fid_vec = &self.group_face_id[&gf];
            for &fid in fid_vec.iter() {
                writeln!(file_group_face, "{},{}", gf, fid).unwrap();
            }
        }

        // return
        Ok(())

    }
}
