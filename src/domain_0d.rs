use crate::domain_1d::Domain1D;
use crate::utils_error::Error1D;

pub struct Domain0D {
    // struct ids
    pub dom0d_id: usize,
    pub dom1d_id: usize,

    // cell data
    pub cell_id: i32,

    // face data
    pub face_id: i32,
    pub face_x: f64,
    pub loc: usize,
}

impl Domain0D {
    pub fn new(dom0d_id: usize, dom1d: &Domain1D, cell_id: i32, loc: usize) -> Result<Domain0D, Error1D> {
        // error checking
        if !dom1d.cell_id.contains(&cell_id) {
            return Err(Error1D::InvalidCellID {caller: "Domain0D".to_string(), cid: cell_id, parent: "Domain1D".to_string()});
        }
        if loc > 1 {
            return Err(Error1D::InvalidLocalID {caller: "Domain0D".to_string(), loc});
        }

        // get struct ids
        let dom1d_id = dom1d.dom1d_id;

        // get face data
        let face_id = dom1d.cell_face_id[&cell_id][loc];
        let face_x = dom1d.face_x[&face_id];

        // return
        Ok(Domain0D {
            dom0d_id,
            dom1d_id,
            cell_id,
            face_id,
            face_x,
            loc,
        })
    }

}