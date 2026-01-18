use crate::problem_1d::Problem1D;
use crate::SteadyBase;
use crate::SteadyDiff;
use crate::SteadyDiffAdv;
use crate::SteadyDiffMulti;
use crate::SteadyDiffAdvMulti;

pub enum SteadyPhysics {
    Diff(SteadyDiff),
    DiffAdv(SteadyDiffAdv),
    DiffMulti(SteadyDiffMulti),
    DiffAdvMulti(SteadyDiffAdvMulti),
}

pub fn parse_steady_physics(
    phys_steady: &mut SteadyPhysics,
    phys_type: &String,
) -> Result<(), crate::utils_error::FVChemError> {

    // initialize physics struct
    match phys_type.as_str() {
        "diff" => {
            *phys_steady = SteadyPhysics::Diff(SteadyDiff::new());
        },
//        "diffadv" => {
//            *phys_steady = Some(Box::new(SteadyDiffAdv::new()));
//        },
//        "diff_multi" => {
//            *phys_steady = Some(Box::new(SteadyDiffMulti::new()));
//        },
//        "diffadv_multi" => {
//            *phys_steady = Some(Box::new(SteadyDiffAdvMulti::new()));
//        },
        _ => {
            panic!("Steady physics type '{}' not recognized.", phys_type);
        },
    }

    // return
    Ok(())

}

pub fn parse_steady_solver(
    prob: &mut Problem1D,
    phys_steady: &mut SteadyPhysics,
    max_iter: usize,
    tol_l2: f64,
    damping: f64,
) -> Result<(), crate::utils_error::FVChemError> {

    // set solver parameters
    match phys_steady {
        SteadyPhysics::Diff(phys) => {
            phys.solve(prob, max_iter, tol_l2, damping)?;
        },
        _ => {
            panic!("Steady physics type not recognized.");
        },
    }

    // return
    Ok(())

}

pub fn parse_steady_condition(
    parts: &Vec<String>,
    phys_steady: &mut SteadyPhysics,
    dom0d_map: &std::collections::HashMap<String, usize>,
    dom1d_map: &std::collections::HashMap<String, usize>,
    scl0d_map: &std::collections::HashMap<String, usize>,
    scl1d_map: &std::collections::HashMap<String, usize>,
    var1d_map: &std::collections::HashMap<String, usize>,
    line_num: usize,
) -> Result<(), crate::utils_error::FVChemError> {

    // TODO: add more physics types 

    // parse dependending on physics type
    match phys_steady {
        SteadyPhysics::Diff(phys) => {
            parse_diff(parts, phys, dom0d_map, dom1d_map, scl0d_map, scl1d_map, var1d_map, line_num)?;
        },
        _ => {
            panic!("Line {}: Steady physics type not recognized.", line_num);
        },
    }

    // return
    Ok(())

}

fn parse_diff(
    parts: &Vec<String>,
    phys: &mut SteadyDiff,
    dom0d_map: &std::collections::HashMap<String, usize>,
    dom1d_map: &std::collections::HashMap<String, usize>,
    scl0d_map: &std::collections::HashMap<String, usize>,
    scl1d_map: &std::collections::HashMap<String, usize>,
    var1d_map: &std::collections::HashMap<String, usize>,
    line_num: usize,
) -> Result<(), crate::utils_error::FVChemError> {

    // TODO: add more physics conditions

    // get inputs
    let cond_type = parts.get(0).expect(format!("Line {}: Expected 'int', 'bnd', or 'itr' for internal, boundary, or interface condition.", line_num).as_str());

    // check if internal, boundary, or interface
    match cond_type.as_str() {
        "int" => {
            // get names of properties
            let dom_name = parts.get(1).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <src>'. Could not find dom.", line_num).as_str());
            let c_name = parts.get(2).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <src>'. Could not find conc.", line_num).as_str());
            let d_name = parts.get(3).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <src>'. Could not find diff.", line_num).as_str());
            let r_name = parts.get(4).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <src>'. Could not find src.", line_num).as_str());

            // get ids
            let dom_id = *dom1d_map.get(dom_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_name).as_str());
            let c_id = *var1d_map.get(c_name).expect(format!("Line {}: var1d name '{}' not found.", line_num, c_name).as_str());
            let d_id = *scl1d_map.get(d_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, d_name).as_str());
            let r_id = *scl1d_map.get(r_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, r_name).as_str());
        
            // add internal domain to physics
            phys.add_domain(dom_id, c_id, d_id, r_id);
        },
        "bnd" => {

            // get boundary type
            let bnd_type = parts.get(1).expect(format!("Line {}: bnd requires 'bnd <bnd_type> ...'. Could not find bnd_type.", line_num).as_str());
            
            // check boundary type
            match bnd_type.as_str() {
                "conc" => {
                    // get names of properties
                    let dom_name = parts.get(2).expect(format!("Line {}: dirichlet requires 'bnd dirichlet <dom> <conc> <value>'. Could not find dom.", line_num).as_str());
                    let c_name = parts.get(3).expect(format!("Line {}: dirichlet requires 'bnd dirichlet <dom> <conc> <value>'. Could not find conc.", line_num).as_str());

                    // get ids
                    let dom_id = *dom0d_map.get(dom_name).expect(format!("Line {}: dom0d name '{}' not found.", line_num, dom_name).as_str());
                    let c_id = *scl0d_map.get(c_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, c_name).as_str());

                    // add dirichlet boundary condition to physics
                    phys.add_boundary_concentration(dom_id, c_id);
                },
                "flux" => {
                    // get names of properties
                    let dom_name = parts.get(2).expect(format!("Line {}: neumann requires 'bnd neumann <dom> <flux> <value>'. Could not find dom.", line_num).as_str());
                    let n_name = parts.get(3).expect(format!("Line {}: neumann requires 'bnd neumann <dom> <flux> <value>'. Could not find flux.", line_num).as_str());

                    // get ids
                    let dom_id = *dom0d_map.get(dom_name).expect(format!("Line {}: dom0d name '{}' not found.", line_num, dom_name).as_str());
                    let n_id = *scl0d_map.get(n_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, n_name).as_str());

                    // add neumann boundary condition to physics
                    phys.add_boundary_flux(dom_id, n_id);
                },
                _ => {
                    panic!("Line {}: Boundary type '{}' not recognized for diffusion physics.", line_num, bnd_type);
                },
            }

        },
        "itr" => {

        },
        _ => {
            panic!("Line {}: Expected 'int', 'bnd', or 'itr' for internal, boundary, or interface condition.", line_num);
        },
    }

    // return
    Ok(())

}
