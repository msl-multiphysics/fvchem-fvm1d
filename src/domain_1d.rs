use crate::mesh_1d::Mesh1D;
use std::collections::{HashMap, HashSet};

pub struct Domain1D {
    // struct ids
    pub dom1d_id: usize,

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
}

impl Domain1D {
    pub fn new(dom1d_id: usize, mesh: &Mesh1D) -> Domain1D {
        // cell data
        let num_cell = mesh.num_cell;
        let cell_id = mesh.cell_id.clone();
        let cell_x = mesh.cell_x.clone();
        let cell_dx = mesh.cell_dx.clone();

        // face data
        let num_face = mesh.num_face;
        let face_id = mesh.face_id.clone();
        let face_x = mesh.face_x.clone();

        // cell to cell data
        let cell_cell_id = mesh.cell_cell_id.clone();
        let cell_cell_dist = mesh.cell_cell_dist.clone();

        // cell to face data
        let cell_face_id = mesh.cell_face_id.clone();
        let cell_face_dist = mesh.cell_face_dist.clone();
        let cell_face_norm = mesh.cell_face_norm.clone();

        // face to cell data
        let face_cell_id = mesh.face_cell_id.clone();
        let face_cell_dist = mesh.face_cell_dist.clone();

        // return
        Domain1D {
            dom1d_id,
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
        }
    }

    pub fn new_from_subset(dom1d_id: usize, mesh: &Mesh1D, cell_id: Vec<i32>) -> Domain1D {
        // cell data
        let num_cell = cell_id.len();
        let mut cell_x: HashMap<i32, f64> = HashMap::new();
        let mut cell_dx: HashMap<i32, f64> = HashMap::new();

        // compute cell data
        for &cid in cell_id.iter() {
            cell_x.insert(cid, mesh.cell_x[&cid]);
            cell_dx.insert(cid, mesh.cell_dx[&cid]);
        }

        // compute face data
        let mut face_id_set: HashSet<i32> = HashSet::new();
        for &cid in cell_id.iter() {
            for loc in 0..2 {
                let fid = mesh.cell_face_id[&cid][loc];
                face_id_set.insert(fid);
            }
        }
        let face_id: Vec<i32> = face_id_set.into_iter().collect();
        let num_face = face_id.len();
        let mut face_x: HashMap<i32, f64> = HashMap::new();
        for &fid in face_id.iter() {
            face_x.insert(fid, mesh.face_x[&fid]);
        }

        // cell to cell data
        let mut cell_cell_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_cell_dist: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to cell data
        for &cid in cell_id.iter() {
            // preliminary neighbor id and distances
            let mut id_sub = mesh.cell_cell_id[&cid].clone();
            let mut dist_sub = mesh.cell_cell_dist[&cid].clone();

            // verify that neighbor cells are in subset
            for loc in 0..2 {
                let nid = mesh.cell_cell_id[&cid][loc]; // neighbor cell id
                if !cell_id.contains(&nid) {
                    // neighbor not in subset
                    id_sub[loc] = mesh.cell_face_id[&cid][loc]; // use face id instead
                    dist_sub[loc] = mesh.cell_face_dist[&cid][loc]; // use face distance instead
                }
            }

            // insert into cell to cell data
            cell_cell_id.insert(cid, id_sub.clone());
            cell_cell_dist.insert(cid, dist_sub.clone());
        }

        // cell to face data
        let mut cell_face_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_face_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        let mut cell_face_norm: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to face data
        for &cid in cell_id.iter() {
            cell_face_id.insert(cid, mesh.cell_face_id[&cid].clone());
            cell_face_dist.insert(cid, mesh.cell_face_dist[&cid].clone());
            cell_face_norm.insert(cid, mesh.cell_face_norm[&cid].clone());
        }

        // face to cell data
        let mut face_cell_id: HashMap<i32, (i32, i32)> = HashMap::new();
        let mut face_cell_dist: HashMap<i32, (f64, f64)> = HashMap::new();

        // compute face to cell data
        for &fid in face_id.iter() {
            // preliminary neighbor ids and distances
            let mut id_sub = mesh.face_cell_id[&fid].clone();
            let mut dist_sub = mesh.face_cell_dist[&fid].clone();

            // verify that neighbor cells are in subset (A)
            let cid_a = mesh.face_cell_id[&fid].0; // neighbor cell A id
            if !cell_id.contains(&cid_a) {
                // neighbor not in subset
                id_sub.0 = fid;
                dist_sub.0 = 0.0;
            }

            // verify that neighbor cells are in subset (B)
            let cid_b = mesh.face_cell_id[&fid].1; // neighbor cell B id
            if !cell_id.contains(&cid_b) {
                // neighbor not in subset
                id_sub.1 = fid;
                dist_sub.1 = 0.0;
            }

            // insert into face to cell data
            face_cell_id.insert(fid, id_sub.clone());
            face_cell_dist.insert(fid, dist_sub.clone());
        }

        // return
        Domain1D {
            dom1d_id,
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
        }
    }
}
