use crate::domain_0d::Domain0D;
use crate::domain_1d::Domain1D;
use crate::mesh_1d::Mesh1D;
use crate::scalar_0d::Scalar0D;
use crate::scalar_1d::Scalar1D;
use crate::variable_1d::Variable1D;

pub struct Problem1D {
    pub dom0d: Vec<Domain0D>,
    pub dom1d: Vec<Domain1D>,
    pub scl0d: Vec<Scalar0D>,
    pub scl1d: Vec<Scalar1D>,
    pub var1d: Vec<Variable1D>,
}

impl Problem1D {
    pub fn new() -> Problem1D {
        Problem1D {
            dom0d: Vec::new(),
            dom1d: Vec::new(),
            scl0d: Vec::new(),
            scl1d: Vec::new(),
            var1d: Vec::new(),
        }
    }

    pub fn add_dom0d(prob: &mut Problem1D, dom1d_id: usize, cell_id: i32, loc: usize) -> usize {
        // get dom0d_id
        let dom0d_id = prob.dom0d.len();

        // create Domain0D
        let dom1d = &prob.dom1d[dom1d_id];
        let dom0d = Domain0D::new(dom0d_id, dom1d, cell_id, loc);
        prob.dom0d.push(dom0d);

        // return
        dom0d_id
    }

    pub fn add_dom1d(prob: &mut Problem1D, mesh: &Mesh1D) -> usize {
        // get dom1d_id
        let dom1d_id = prob.dom1d.len();

        // create Domain1D
        let dom1d = Domain1D::new(dom1d_id, mesh);
        prob.dom1d.push(dom1d);

        // return
        dom1d_id
    }

    pub fn add_dom1d_from_subset(prob: &mut Problem1D, mesh: &Mesh1D, cell_id: Vec<i32>) -> usize {
        // get dom1d_id
        let dom1d_id = prob.dom1d.len();

        // create Domain1D
        let dom1d = Domain1D::new_from_subset(dom1d_id, mesh, cell_id);
        prob.dom1d.push(dom1d);

        // return
        dom1d_id
    }

    pub fn add_scl0d(
        prob: &mut Problem1D,
        dom0d_id: usize,
        value: f64,
        output_file: String,
        output_step: usize,
    ) -> usize {
        // get scl0d_id
        let scl0d_id = prob.scl0d.len();

        // create Scalar0D
        let dom0d = &prob.dom0d[dom0d_id];
        let scl0d = Scalar0D::new(scl0d_id, dom0d, value, output_file, output_step);
        prob.scl0d.push(scl0d);

        // return
        scl0d_id
    }

    pub fn add_scl0d_nonconstant(
        prob: &mut Problem1D,
        dom0d_id: usize,
        value_func: fn(f64, f64, Vec<f64>) -> f64,
        var1d_id: Vec<usize>,
        output_file: String,
        output_step: usize,
    ) -> usize {
        // get scl0d_id
        let scl0d_id = prob.scl0d.len();

        // create Scalar0D
        let dom0d = &prob.dom0d[dom0d_id];
        let scl0d = Scalar0D::new_nonconstant(
            scl0d_id,
            dom0d,
            &prob.var1d,
            var1d_id,
            value_func,
            output_file,
            output_step,
        );
        prob.scl0d.push(scl0d);

        // return
        scl0d_id
    }

    pub fn add_scl1d(
        prob: &mut Problem1D,
        dom1d_id: usize,
        value: f64,
        output_file: String,
        output_step: usize,
    ) -> usize {
        // get scl1d_id
        let scl1d_id = prob.scl1d.len();

        // create Scalar1D
        let dom1d = &prob.dom1d[dom1d_id];
        let scl1d = Scalar1D::new(scl1d_id, dom1d, value, output_file, output_step);
        prob.scl1d.push(scl1d);

        // return
        scl1d_id
    }

    pub fn add_scl1d_nonconstant(
        prob: &mut Problem1D,
        dom1d_id: usize,
        value_func: fn(f64, f64, Vec<f64>) -> f64,
        var1d_id: Vec<usize>,
        output_file: String,
        output_step: usize,
    ) -> usize {
        // get scl1d_id
        let scl1d_id = prob.scl1d.len();

        // create Scalar1D
        let dom1d = &prob.dom1d[dom1d_id];
        let scl1d = Scalar1D::new_nonconstant(
            scl1d_id,
            dom1d,
            &prob.var1d,
            var1d_id,
            value_func,
            output_file,
            output_step,
        );
        prob.scl1d.push(scl1d);

        // return
        scl1d_id
    }

    pub fn add_var1d(
        prob: &mut Problem1D,
        dom1d_id: usize,
        value_init: f64,
        output_file: String,
        output_step: usize,
    ) -> usize {
        // get var1d_id
        let var1d_id = prob.var1d.len();

        // create Variable1D
        let dom1d = &prob.dom1d[dom1d_id];
        let var1d = Variable1D::new(var1d_id, dom1d, value_init, output_file, output_step);
        prob.var1d.push(var1d);

        // return
        var1d_id
    }
}
