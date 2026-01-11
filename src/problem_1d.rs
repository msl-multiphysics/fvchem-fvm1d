use crate::domain_0d::Domain0D;
use crate::domain_1d::Domain1D;
use crate::utils_error::Error1D;
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

    pub fn add_dom0d(&mut self, dom1d_id: usize, cell_id: i32, loc: usize) -> Result<usize, Error1D> {
        // get dom0d_id
        let dom0d_id = self.dom0d.len();

        // create Domain0D
        let dom1d = &self.dom1d[dom1d_id];
        let dom0d = Domain0D::new(dom0d_id, dom1d, cell_id, loc)?;
        self.dom0d.push(dom0d);

        // return
        Ok(dom0d_id)
    }

    pub fn add_dom1d(&mut self, mesh: &Mesh1D) -> Result<usize, Error1D> {
        // get dom1d_id
        let dom1d_id = self.dom1d.len();

        // create Domain1D
        let dom1d = Domain1D::new(dom1d_id, mesh)?;
        self.dom1d.push(dom1d);

        // return
        Ok(dom1d_id)
    }

    pub fn add_dom1d_from_subset(&mut self, mesh: &Mesh1D, cell_id: Vec<i32>) -> Result<usize, Error1D> {
        // get dom1d_id
        let dom1d_id = self.dom1d.len();

        // create Domain1D
        let dom1d = Domain1D::new_from_subset(dom1d_id, mesh, cell_id)?;
        self.dom1d.push(dom1d);

        // return
        Ok(dom1d_id)
    }

    pub fn add_scl0d(
        &mut self,
        dom0d_id: usize,
        value: f64
    ) -> Result<usize, Error1D> {
        // get scl0d_id
        let scl0d_id = self.scl0d.len();

        // create Scalar0D
        let dom0d = &self.dom0d[dom0d_id];
        let scl0d = Scalar0D::new(scl0d_id, dom0d, value)?;
        self.scl0d.push(scl0d);

        // return
        Ok(scl0d_id)
    }

    pub fn add_scl0d_from_function(
        &mut self,
        dom0d_id: usize,
        value_func: fn(usize, f64, Vec<f64>) -> f64,
        value_var: Vec<usize>,
    ) -> Result<usize, Error1D> {
        // get scl0d_id
        let scl0d_id = self.scl0d.len();

        // create Scalar0D
        let dom0d = &self.dom0d[dom0d_id];
        let scl0d = Scalar0D::new_from_function(
            scl0d_id,
            dom0d,
            &self.var1d,
            value_func,
            value_var,
        )?;
        self.scl0d.push(scl0d);

        // return
        Ok(scl0d_id)
    }

    pub fn add_scl0d_from_file(
        &mut self,
        dom0d_id: usize,
        value_file: String,
    ) -> Result<usize, Error1D> {
        // get scl0d_id
        let scl0d_id = self.scl0d.len();

        // create Scalar0D
        let dom0d = &self.dom0d[dom0d_id];
        let scl0d = Scalar0D::new_from_file(scl0d_id, dom0d, value_file)?;
        self.scl0d.push(scl0d);

        // return
        Ok(scl0d_id)
    }

    pub fn add_scl1d(
        &mut self,
        dom1d_id: usize,
        value: f64,
    ) -> Result<usize, Error1D> {
        // get scl1d_id
        let scl1d_id = self.scl1d.len();

        // create Scalar1D
        let dom1d = &self.dom1d[dom1d_id];
        let scl1d = Scalar1D::new(scl1d_id, dom1d, value)?;
        self.scl1d.push(scl1d);

        // return
        Ok(scl1d_id)
    }

    pub fn add_scl1d_from_function(
        &mut self,
        dom1d_id: usize,
        value_func: fn(usize, f64, Vec<f64>) -> f64,
        value_var: Vec<usize>,
    ) -> Result<usize, Error1D> {
        // get scl1d_id
        let scl1d_id = self.scl1d.len();

        // create Scalar1D
        let dom1d = &self.dom1d[dom1d_id];
        let scl1d = Scalar1D::new_from_function(
            scl1d_id,
            dom1d,
            &self.var1d,
            value_func,
            value_var,
        )?;
        self.scl1d.push(scl1d);
        
        // return
        Ok(scl1d_id)
    }

    pub fn add_scl1d_from_file(
        &mut self,
        dom1d_id: usize,
        value_file: String,
    ) -> Result<usize, Error1D> {
        // get scl1d_id
        let scl1d_id = self.scl1d.len();
        
        // create Scalar1D
        let dom1d = &self.dom1d[dom1d_id];
        let scl1d = Scalar1D::new_from_file(scl1d_id, dom1d, value_file)?;
        self.scl1d.push(scl1d);

        // return
        Ok(scl1d_id)
    }

    pub fn add_var1d(
        &mut self,
        dom1d_id: usize,
        value_init: f64,
    ) -> Result<usize, Error1D> {
        // get var1d_id
        let var1d_id = self.var1d.len();

        // create Variable1D
        let dom1d = &self.dom1d[dom1d_id];
        let var1d = Variable1D::new(var1d_id, dom1d, value_init)?;
        self.var1d.push(var1d);

        // return
        Ok(var1d_id)
    }

    pub fn set_scl0d_write_steady(&mut self, scl0d_id: usize, write_file: String) {
        Scalar0D::set_write_steady(&mut self.scl0d[scl0d_id], write_file);
    }

    pub fn set_scl0d_write_transient(&mut self, scl0d_id: usize, write_file: String, write_step: usize) {
        Scalar0D::set_write_transient(&mut self.scl0d[scl0d_id], write_file, write_step);
    }

    pub fn set_scl1d_write_steady(&mut self, scl1d_id: usize, write_file: String) {
        Scalar1D::set_write_steady(&mut self.scl1d[scl1d_id], write_file);
    }

    pub fn set_scl1d_write_transient(&mut self, scl1d_id: usize, write_file: String, write_step: usize) {
        Scalar1D::set_write_transient(&mut self.scl1d[scl1d_id], write_file, write_step);
    }

    pub fn set_var1d_write_steady(&mut self, var1d_id: usize, write_file: String) {
        Variable1D::set_write_steady(&mut self.var1d[var1d_id], write_file);
    }

    pub fn set_var1d_write_transient(&mut self, var1d_id: usize, write_file: String, write_step: usize) {
        Variable1D::set_write_transient(&mut self.var1d[var1d_id], write_file, write_step);
    }

}
