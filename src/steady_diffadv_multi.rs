use crate::problem_1d::Problem1D;
use crate::steady_base::SteadyBase;
use crate::utils_limiter::{LimiterType, calc_limiter};
use faer::prelude::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

pub struct SteadyDiffAdvMulti {
    // internal data
    // key: (domain id, component i id)
    pub internal_dom: Vec<(usize, usize)>,
    pub internal_c: HashMap<(usize, usize), usize>,  // c_i
    pub internal_d: HashMap<(usize, usize), HashMap<usize, usize>>,  // D_ij wherein inner key is component j
    pub internal_u: HashMap<(usize, usize), usize>,  // u_i
    pub internal_r: HashMap<(usize, usize), usize>,  // R_i

    // boundary data
    // key: (domain id, component id)
    pub bndconc_dom: Vec<(usize, usize)>, // concentration boundary
    pub bndconc_c: HashMap<(usize, usize), usize>,
    pub bndflux_dom: Vec<(usize, usize)>, // flux boundary
    pub bndflux_n: HashMap<(usize, usize), usize>,
    pub bndnadv_dom: Vec<(usize, usize)>, // flux (not counting advection) boundary
    pub bndnadv_n: HashMap<(usize, usize), usize>,
    pub bndmtrn_dom: Vec<(usize, usize)>, // mass transfer boundary
    pub bndmtrn_k: HashMap<(usize, usize), usize>,
    pub bndmtrn_c: HashMap<(usize, usize), usize>,

    // interface data
    // key: (domain a, domain b, component id)
    pub itrcont_dom: Vec<(usize, usize, usize)>, // continuity interface
    pub itrmtrn_dom: Vec<(usize, usize, usize)>, // mass transfer interface
    pub itrmtrn_k: HashMap<(usize, usize, usize), usize>,

    // number of components and limiter type
    pub num_comp: usize,
    pub limiter_type: LimiterType,
}

impl SteadyBase for SteadyDiffAdvMulti {
    fn assemble_matrix(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
    ) {
        // loop over internal domains
        for &(dom1d_id, comp_id) in self.internal_dom.iter() {
            self.assemble_internal(prob, a_triplet, b_vec, dom1d_id, comp_id);
        }

        // loop over boundary conditions
        for &(dom0d_id, comp_id) in self.bndconc_dom.iter() {
            self.assemble_bndconc(prob, a_triplet, b_vec, dom0d_id, comp_id);
        }
        for &(dom0d_id, comp_id) in self.bndflux_dom.iter() {
            self.assemble_bndflux(prob, a_triplet, b_vec, dom0d_id, comp_id);
        }
        for &(dom0d_id, comp_id) in self.bndnadv_dom.iter() {
            self.assemble_bndnadv(prob, a_triplet, b_vec, dom0d_id, comp_id);
        }   
        for &(dom0d_id, comp_id) in self.bndmtrn_dom.iter() {
            self.assemble_bndmtrn(prob, a_triplet, b_vec, dom0d_id, comp_id);
        }

        // loop over interface conditions
        for &(dom0d_a, dom0d_b, comp_id) in self.itrcont_dom.iter() {
            self.assemble_itrcont(prob, a_triplet, dom0d_a, dom0d_b, comp_id);
        }
        for &(dom0d_a, dom0d_b, comp_id) in self.itrmtrn_dom.iter() {
            self.assemble_itrmtrn(prob, a_triplet, dom0d_a, dom0d_b, comp_id);
        }
    }
}

impl SteadyDiffAdvMulti {
    pub fn new(num_comp: usize, limiter_type: LimiterType) -> SteadyDiffAdvMulti {
        SteadyDiffAdvMulti {
            internal_dom: Vec::new(),
            internal_c: HashMap::new(),
            internal_d: HashMap::new(),
            internal_u: HashMap::new(),
            internal_r: HashMap::new(),
            bndconc_dom: Vec::new(),
            bndconc_c: HashMap::new(),
            bndflux_dom: Vec::new(),
            bndflux_n: HashMap::new(),
            bndnadv_dom: Vec::new(),
            bndnadv_n: HashMap::new(),
            bndmtrn_dom: Vec::new(),
            bndmtrn_k: HashMap::new(),
            bndmtrn_c: HashMap::new(),
            itrcont_dom: Vec::new(),
            itrmtrn_dom: Vec::new(),
            itrmtrn_k: HashMap::new(),
            num_comp,
            limiter_type,
        }
    }

    pub fn add_domain(&mut self, dom1d: usize, comp_id: usize, conc: usize, diff_coeff: HashMap<usize, usize>, src: usize) {
        self.internal_dom.push((dom1d, comp_id));
        self.internal_c.insert((dom1d, comp_id), conc);
        self.internal_d.insert((dom1d, comp_id), diff_coeff);
        self.internal_r.insert((dom1d, comp_id), src);
    }

    pub fn add_boundary_concentration(&mut self, dom0d: usize, comp_id: usize, conc: usize) {
        self.bndconc_dom.push((dom0d, comp_id));
        self.bndconc_c.insert((dom0d, comp_id), conc);
    }

    pub fn add_boundary_flux(&mut self, dom0d: usize, comp_id: usize, flux: usize) {
        self.bndflux_dom.push((dom0d, comp_id));
        self.bndflux_n.insert((dom0d, comp_id), flux);
    }

    pub fn add_boundary_fluxnoadvection(&mut self, dom0d: usize, comp_id: usize, flux: usize) {
        self.bndnadv_dom.push((dom0d, comp_id));
        self.bndnadv_n.insert((dom0d, comp_id), flux);
    }

    pub fn add_boundary_masstransfer(&mut self, dom0d: usize, comp_id: usize, mass_coeff: usize, conc_ref: usize) {
        self.bndmtrn_dom.push((dom0d, comp_id));
        self.bndmtrn_k.insert((dom0d, comp_id), mass_coeff);
        self.bndmtrn_c.insert((dom0d, comp_id), conc_ref);
    }

    pub fn add_interface_continuity(&mut self, dom0d_a: usize, dom0d_b: usize, comp_id: usize) {
        self.itrcont_dom.push((dom0d_a, dom0d_b, comp_id));
    }

    pub fn add_interface_masstransfer(
        &mut self,
        dom0d_a: usize,
        dom0d_b: usize,
        comp_id: usize,
        mass_coeff: usize,
    ) {
        self.itrmtrn_dom.push((dom0d_a, dom0d_b, comp_id));
        self.itrmtrn_k.insert((dom0d_a, dom0d_b, comp_id), mass_coeff);
    }

    fn assemble_internal(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dom1d_id: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let c_id = self.internal_c[&(dom1d_id, comp_id)];

        // loop over cells
        for &cid in prob.dom1d[dom1d_id].cell_id.iter() {
            // flux term
            for loc in 0..2 {
                let fid = prob.dom1d[dom1d_id].cell_face_id[&cid][loc];
                let nid = prob.dom1d[dom1d_id].cell_cell_id[&cid][loc];
                if nid >= 0 {
                    // internal face
                    self.assemble_flux_cn(prob, a_triplet, dom1d_id, comp_id, true, c_id, cid, cid, fid, nid, loc); // discretized equation - store in cell equation
                    self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, true, c_id, fid, cid, fid, loc); // face interpolation - store in face equation
                    self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, true, c_id, fid, nid, fid, loc); // face interpolation - store in face equation
                } else {
                    self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, true, c_id, cid, cid, fid, loc); // boundary discretized equation - store in cell equation
                    // face equation handled in boundary and interface conditions
                }
            }

            // source term
            self.assemble_src(prob, b_vec, dom1d_id, comp_id, cid, cid); // source term - store in cell equation
        }
    }

    fn assemble_flux_cn(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        dom1d_id: usize,
        comp_id: usize,
        with_adv: bool,
        var_row: usize,
        row: i32,
        cid: i32,
        fid: i32,
        nid: i32,
        loc: usize,
    ) {
        // diffusion term
        
        // this component (comp_id)
        let d_map = &self.internal_d[&(dom1d_id, comp_id)];
        
        // get properties
        let dist_cn = prob.dom1d[dom1d_id].cell_cell_dist[&cid][loc];

        // diffusion coefficients affecting this component (comp_other)
        for (&comp_other, &d_id) in d_map.iter() {
            // get variable ids
            let c_other = self.internal_c[&(dom1d_id, comp_other)];

            // get properties
            let d_f = prob.scl1d[d_id].face_value[&fid];
            
            // add to matrix
            self.add_a(prob, a_triplet, var_row, c_other, row, cid,  d_f / dist_cn);
            self.add_a(prob, a_triplet, var_row, c_other, row, nid, -d_f / dist_cn);
        }

        // advection term

        // skip if no advection
        if !with_adv {
            return;
        }

        // get variable ids
        let c_id = self.internal_c[&(dom1d_id, comp_id)];
        let u_id = self.internal_u[&(dom1d_id, comp_id)];

        // get properties
        let u_f = prob.scl1d[u_id].face_value[&fid];
        let norm_cf = prob.dom1d[dom1d_id].cell_face_norm[&cid][loc];
        let dist_cf = prob.dom1d[dom1d_id].cell_face_dist[&cid][loc];
        let u_norm = u_f * norm_cf;
        let limiter = calc_limiter(prob, dom1d_id, c_id, cid, &self.limiter_type);

        // identify upwind and downwind cells
        let uid = if u_norm >= 0.0 { cid } else { nid }; // upwind
        let did = if u_norm >= 0.0 { nid } else { cid }; // downwind

        // advection term
        self.add_a(
            prob,
            a_triplet,
            var_row,
            c_id,
            row,
            uid,
            u_norm * (1.0 - limiter * dist_cf / dist_cn),
        );
        self.add_a(
            prob,
            a_triplet,
            var_row,
            c_id,
            row,
            did,
            u_norm * limiter * dist_cf / dist_cn,
        );
    }

    fn assemble_flux_cf(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        dom1d_id: usize,
        comp_id: usize,
        with_adv: bool,
        var_row: usize,
        row: i32,
        cid: i32,
        fid: i32,
        loc: usize,
    ) {
        // diffusion term

        // this component (comp_row)
        let d_map = &self.internal_d[&(dom1d_id, comp_id)];

        // get properties
        let dist_cf = prob.dom1d[dom1d_id].cell_face_dist[&cid][loc];

        // diffusion coefficients affecting this component (comp_other)
        for (&comp_other, &d_id) in d_map.iter() {
            // get variable ids
            let c_other = self.internal_c[&(dom1d_id, comp_other)];

            // get properties
            let d_f = prob.scl1d[d_id].face_value[&fid];
            
            // add to matrix
            self.add_a(prob, a_triplet, var_row, c_other, row, cid, d_f / dist_cf);
            self.add_a(prob, a_triplet, var_row, c_other, row, fid, -d_f / dist_cf);
        }

        // advection term

        // skip if no advection
        if !with_adv {
            return;
        }

        // get variable ids
        let c_id = self.internal_c[&(dom1d_id, comp_id)];
        let u_id = self.internal_u[&(dom1d_id, comp_id)];

        // get properties
        let u_f = prob.scl1d[u_id].face_value[&fid];
        let norm_cf = prob.dom1d[dom1d_id].cell_face_norm[&cid][loc];
        let u_norm = u_f * norm_cf;
        let limiter = calc_limiter(prob, dom1d_id, c_id, cid, &self.limiter_type);

        // identify upwind and downwind cells
        let uid = if u_norm >= 0.0 { cid } else { fid }; // upwind
        let did = if u_norm >= 0.0 { fid } else { cid }; // downwind

        // advection term
        self.add_a(
            prob,
            a_triplet,
            var_row,
            c_id,
            row,
            uid,
            u_norm * (1.0 - limiter),
        );
        self.add_a(prob, a_triplet, var_row, c_id, row, did, u_norm * limiter);
    }

    fn assemble_src(
        &self,
        prob: &Problem1D,
        b_vec: &mut Col<f64>,
        dom1d_id: usize,
        comp_id: usize,
        row: i32,
        cid: i32,
    ) {
        // this component (comp_id)
        let c_id = self.internal_c[&(dom1d_id, comp_id)];
        let r_id = self.internal_r[&(dom1d_id, comp_id)];

        // get properties
        let r_c = prob.scl1d[r_id].cell_value[&cid];
        let dx_c = prob.dom1d[dom1d_id].cell_dx[&cid];

        // add to rhs
        self.add_b(prob, b_vec, c_id, row, r_c * dx_c);
    }

    fn assemble_bndconc(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dom0d_id: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_id = prob.dom0d[dom0d_id].dom1d_id;
        let c0d_id = self.bndconc_c[&(dom0d_id, comp_id)];
        let c_id = self.internal_c[&(dom1d_id, comp_id)];

        // get properties
        let fid = prob.dom0d[dom0d_id].face_id;
        let c_f = prob.scl0d[c0d_id].face_value;

        // add to matrix
        self.add_a(prob, a_triplet, c_id, c_id, fid, fid, 1.0);
        self.add_b(prob, b_vec, c_id, fid, c_f);
    }

    fn assemble_bndflux(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dom0d_id: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_id = prob.dom0d[dom0d_id].dom1d_id;
        let n0d_id = self.bndflux_n[&(dom0d_id, comp_id)];
        let c_id = self.internal_c[&(dom1d_id, comp_id)];

        // get properties
        let cid = prob.dom0d[dom0d_id].cell_id;
        let fid = prob.dom0d[dom0d_id].face_id;
        let loc = prob.dom0d[dom0d_id].loc;
        let n_f = prob.scl0d[n0d_id].face_value;

        // add to rhs
        self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, true, c_id, fid, cid, fid, loc); // flux discretization
        self.add_b(prob, b_vec, c_id, fid, n_f);
    }

    fn assemble_bndnadv(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dom0d_id: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_id = prob.dom0d[dom0d_id].dom1d_id;
        let n0d_id = self.bndnadv_n[&(dom0d_id, comp_id)];
        let c_id = self.internal_c[&(dom1d_id, comp_id)];

        // get properties
        let cid = prob.dom0d[dom0d_id].cell_id;
        let fid = prob.dom0d[dom0d_id].face_id;
        let loc = prob.dom0d[dom0d_id].loc;
        let n_f = prob.scl0d[n0d_id].face_value;

        // add to rhs
        self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, false, c_id, fid, cid, fid, loc); // flux discretization
        self.add_b(prob, b_vec, c_id, fid, n_f);
    }

    fn assemble_bndmtrn(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        b_vec: &mut Col<f64>,
        dom0d_id: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_id = prob.dom0d[dom0d_id].dom1d_id;
        let k0d_id = self.bndmtrn_k[&(dom0d_id, comp_id)];
        let c0d_id = self.bndmtrn_c[&(dom0d_id, comp_id)];
        let c_id = self.internal_c[&(dom1d_id, comp_id)];

        // get properties
        let cid = prob.dom0d[dom0d_id].cell_id;
        let fid = prob.dom0d[dom0d_id].face_id;
        let loc = prob.dom0d[dom0d_id].loc;
        let k_f = prob.scl0d[k0d_id].face_value;
        let c_ref = prob.scl0d[c0d_id].face_value;

        // add to matrix and rhs
        self.assemble_flux_cf(prob, a_triplet, dom1d_id, comp_id, false, c_id, fid, cid, fid, loc); // flux discretization
        self.add_a(prob, a_triplet, c_id, c_id, fid, fid, -k_f);
        self.add_b(prob, b_vec, c_id, fid, -k_f * c_ref);
    }

    fn assemble_itrcont(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        dom0d_a: usize,
        dom0d_b: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_a = prob.dom0d[dom0d_a].dom1d_id;
        let dom1d_b = prob.dom0d[dom0d_b].dom1d_id;
        let c_id_a = self.internal_c[&(dom1d_a, comp_id)];
        let c_id_b = self.internal_c[&(dom1d_b, comp_id)];

        // get properties
        let cid_a = prob.dom0d[dom0d_a].cell_id;
        let cid_b = prob.dom0d[dom0d_b].cell_id;
        let fid_a = prob.dom0d[dom0d_a].face_id;
        let fid_b = prob.dom0d[dom0d_b].face_id;
        let loc_a = prob.dom0d[dom0d_a].loc;
        let loc_b = prob.dom0d[dom0d_b].loc;

        // concentration continuity - store in face A
        self.add_a(prob, a_triplet, c_id_a, c_id_a, fid_a, fid_a, 1.0);
        self.add_a(prob, a_triplet, c_id_a, c_id_b, fid_a, fid_b, -1.0);

        // flux continuity - store in face B
        self.assemble_flux_cf(prob, a_triplet, dom1d_a, comp_id, true, c_id_b, fid_b, cid_a, fid_a, loc_a);
        self.assemble_flux_cf(prob, a_triplet, dom1d_b, comp_id, true, c_id_b, fid_b, cid_b, fid_b, loc_b);
    }

    fn assemble_itrmtrn(
        &self,
        prob: &Problem1D,
        a_triplet: &mut Vec<Triplet<usize, usize, f64>>,
        dom0d_a: usize,
        dom0d_b: usize,
        comp_id: usize,
    ) {
        // get variable ids
        let dom1d_a = prob.dom0d[dom0d_a].dom1d_id;
        let dom1d_b = prob.dom0d[dom0d_b].dom1d_id;
        let c_id_a = self.internal_c[&(dom1d_a, comp_id)];
        let c_id_b = self.internal_c[&(dom1d_b, comp_id)];
        let k_id = self.itrmtrn_k[&(dom0d_a, dom0d_b, comp_id)];

        // get properties
        let cid_a = prob.dom0d[dom0d_a].cell_id;
        let cid_b = prob.dom0d[dom0d_b].cell_id;
        let fid_a = prob.dom0d[dom0d_a].face_id;
        let fid_b = prob.dom0d[dom0d_b].face_id;
        let loc_a = prob.dom0d[dom0d_a].loc;
        let loc_b = prob.dom0d[dom0d_b].loc;
        let k_f = prob.scl0d[k_id].face_value;

        // concentration jump - store in face A
        self.assemble_flux_cf(prob, a_triplet, dom1d_a, comp_id, false, c_id_a, fid_a, cid_a, fid_a, loc_a);
        self.add_a(prob, a_triplet, c_id_a, c_id_a, fid_a, fid_a, -k_f);
        self.add_a(prob, a_triplet, c_id_a, c_id_b, fid_a, fid_b, k_f);

        // flux continuity - store in face B
        self.assemble_flux_cf(prob, a_triplet, dom1d_a, comp_id, true, c_id_b, fid_b, cid_a, fid_a, loc_a);
        self.assemble_flux_cf(prob, a_triplet, dom1d_b, comp_id, true, c_id_b, fid_b, cid_b, fid_b, loc_b);
    }
}
