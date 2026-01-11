use crate::error_1d::Error1D;
use std::collections::HashMap;

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
        })
    }

    pub fn new_from_face(face_x_vec: Vec<f64>) -> Result<Mesh1D, Error1D> {
        // error checking
        if face_x_vec.len() < 2 {
            return Err(Error1D::InvalidCellCount {caller: "Mesh1D".to_string()});
        }

        // sort face positions
        let mut face_x_sort = face_x_vec.clone();
        face_x_sort.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // face data
        let mut face_id: Vec<i32> = Vec::new();
        let mut face_x: HashMap<i32, f64> = HashMap::new();

        // compute face data
        let num_face = face_x_sort.len();
        for i in 0..num_face {
            let fid = -(i as i32 + 1); // face id: -1 to -num_face
            face_id.push(fid);
            face_x.insert(fid, face_x_sort[i]);
        }

        // cell data
        let mut cell_id: Vec<i32> = Vec::new();
        let mut cell_x: HashMap<i32, f64> = HashMap::new();
        let mut cell_dx: HashMap<i32, f64> = HashMap::new();

        // compute cell data
        let num_cell = num_face - 1;
        for i in 0..num_cell {
            let cid = i as i32; // cell id: 0 to num_cell-1
            let dx = face_x_sort[i + 1] - face_x_sort[i];
            cell_id.push(cid);
            cell_x.insert(cid, 0.5 * (face_x_sort[i] + face_x_sort[i + 1]));
            cell_dx.insert(cid, dx);
            if dx <= 0.0 {
                return Err(Error1D::InvalidCellSize {caller: "Mesh1D".to_string()});
            }
        }

        // cell to cell data
        let mut cell_cell_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_cell_dist: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to cell data
        for i in 1..num_cell - 1 {
            // interior cells
            let cid = i as i32;
            cell_cell_id.insert(cid, vec![i as i32 - 1, i as i32 + 1]);
            cell_cell_dist.insert(
                cid,
                vec![
                    cell_x[&(i as i32)] - cell_x[&(i as i32 - 1)],
                    cell_x[&(i as i32 + 1)] - cell_x[&(i as i32)],
                ],
            );
        }
        cell_cell_id.insert(0, vec![-1, 1]); // no neighbor cell -> default to face
        cell_cell_dist.insert(0, vec![face_x[&(-1)] - cell_x[&0], cell_x[&1] - cell_x[&0]]);
        cell_cell_id.insert(
            num_cell as i32 - 1,
            vec![(num_cell as i32 - 2), -(num_face as i32)],
        );
        cell_cell_dist.insert(
            num_cell as i32 - 1,
            vec![
                cell_x[&(num_cell as i32 - 1)] - cell_x[&(num_cell as i32 - 2)],
                face_x[&(-(num_face as i32))] - cell_x[&(num_cell as i32 - 1)],
            ],
        );

        // cell to face data
        let mut cell_face_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut cell_face_dist: HashMap<i32, Vec<f64>> = HashMap::new();
        let mut cell_face_norm: HashMap<i32, Vec<f64>> = HashMap::new();

        // compute cell to face data
        for i in 0..num_cell {
            let cid = i as i32;
            cell_face_id.insert(cid, vec![-(i as i32 + 1), -(i as i32 + 2)]);
            cell_face_dist.insert(
                cid,
                vec![
                    cell_x[&cid] - face_x[&(-(i as i32 + 1))],
                    face_x[&(-(i as i32 + 2))] - cell_x[&cid],
                ],
            );
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
            face_cell_dist.insert(
                fid,
                (
                    face_x[&fid] - cell_x[&(i as i32 - 1)],
                    cell_x[&(i as i32)] - face_x[&fid],
                ),
            );
        }
        face_cell_id.insert(-1, (-1, 0)); // no neighbor cell -> default to face
        face_cell_dist.insert(-1, (0.0, cell_x[&0] - face_x[&(-1)]));
        face_cell_id.insert(
            -(num_face as i32),
            ((num_cell as i32 - 1), -(num_face as i32)),
        );
        face_cell_dist.insert(
            -(num_face as i32),
            (
                face_x[&(-(num_face as i32))] - cell_x[&(num_cell as i32 - 1)],
                0.0,
            ),
        );

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
        })
    }

}
