use crate::problem_1d::Problem1D;
use crate::utils_error::FVChemError;
use crate::utils_limiter::LimiterType;

use crate::TransientBase;
use crate::TransientDiff;
use crate::TransientDiffAdv;
use crate::TransientDiffMulti;
use crate::TransientDiffAdvMulti;

pub enum TransientPhysics {
    Diff(TransientDiff),
    DiffAdv(TransientDiffAdv),
    DiffMulti(TransientDiffMulti),
    DiffAdvMulti(TransientDiffAdvMulti),
}

pub fn parse_transient_physics(
    parts: &Vec<String>,
    phys_transient: &mut Option<TransientPhysics>,
    phys_type: &String,
    line_num: usize,
) -> Result<(), FVChemError> {

    // initialize physics struct
    match phys_type.as_str() {
        "diff" => {
            *phys_transient = Some(TransientPhysics::Diff(TransientDiff::new()));
        },
        "diff_multi" => {
            let num_comp = parts.get(2).expect(format!("Line {}: physics diff_multi requires 'physics diff_multi <num_comp>'. Could not find num_comp.", line_num).as_str()).parse().expect(format!("Line {}: Could not parse num_comp into a number.", line_num).as_str());
            *phys_transient = Some(TransientPhysics::DiffMulti(TransientDiffMulti::new(num_comp)));
        },
        "diffadv" => {
            // get input
            let limit_str = parts.get(2).expect(format!("Line {}: physics diffadv requires 'physics diffadv <limiter>'. Could not find limiter.", line_num).as_str());
            
            // get limiter type
            let limit_type = match limit_str.as_str() {
                "linear" | "central" | "lin" => LimiterType::Linear,
                "upwind" => LimiterType::Upwind,
                "barth_jespersen" | "barth-jespersen" | "barthjespersen" | "bj" => LimiterType::BarthJespersen,
                "venkatakrishnan" | "ven" | "vk" => LimiterType::Venkatakrishnan,
                _ => {
                    panic!("Line {}: Limiter type '{}' not recognized.", line_num, limit_str);
                },
            };
            
            // set physics
            *phys_transient = Some(TransientPhysics::DiffAdv(TransientDiffAdv::new(limit_type)));
        },
        "diffadv_multi" => {
            // get input
            let num_comp = parts.get(2).expect(format!("Line {}: physics diffadv_multi requires 'physics diffadv_multi <limiter> <num_comp>'. Could not find num_comp.", line_num).as_str()).parse().expect(format!("Line {}: Could not parse num_comp into a number.", line_num).as_str());
            let limit_str = parts.get(3).expect(format!("Line {}: physics diffadv_multi requires 'physics diffadv_multi <limiter>'. Could not find limiter.", line_num).as_str());
            
            // get limiter type
            let limit_type = match limit_str.as_str() {
                "linear" | "central" | "lin" => LimiterType::Linear,
                "upwind" => LimiterType::Upwind,
                "barth_jespersen" | "barth-jespersen" | "barthjespersen" | "bj" => LimiterType::BarthJespersen,
                "venkatakrishnan" | "ven" | "vk" => LimiterType::Venkatakrishnan,
                _ => {
                    panic!("Line {}: Limiter type '{}' not recognized.", line_num, limit_str);
                },
            };

            // set physics
            *phys_transient = Some(TransientPhysics::DiffAdvMulti(TransientDiffAdvMulti::new(num_comp, limit_type)));
        },
        _ => {
            panic!("Transient physics type '{}' not recognized.", phys_type);
        },
    }

    // return
    Ok(())

}

pub fn parse_transient_solver(
    prob: &mut Problem1D,
    phys_transient: &mut Option<TransientPhysics>,
    dt: f64,
    num_step: usize,
    max_iter: usize,
    tol_l2: f64,
    damping: f64,
) -> Result<(), FVChemError> {

    // set solver parameters
    match phys_transient {
        Some(TransientPhysics::Diff(phys)) => {
            phys.solve(prob, dt, num_step, max_iter, tol_l2, damping)?;
        },
        Some(TransientPhysics::DiffAdv(phys)) => {
            phys.solve(prob, dt, num_step, max_iter, tol_l2, damping)?;
        },
        Some(TransientPhysics::DiffMulti(phys)) => {
            phys.solve(prob, dt, num_step, max_iter, tol_l2, damping)?;
        },
        Some(TransientPhysics::DiffAdvMulti(phys)) => {
            phys.solve(prob, dt, num_step, max_iter, tol_l2, damping)?;
        },
        _ => {
            panic!("Transient physics type not recognized.");
        },
    }

    // return
    Ok(())

}

pub fn parse_transient_condition(
    parts: &Vec<String>,
    phys_transient: &mut Option<TransientPhysics>,
    dom0d_map: &std::collections::HashMap<String, usize>,
    dom1d_map: &std::collections::HashMap<String, usize>,
    scl0d_map: &std::collections::HashMap<String, usize>,
    scl1d_map: &std::collections::HashMap<String, usize>,
    var1d_map: &std::collections::HashMap<String, usize>,
    line_num: usize,
) -> Result<(), FVChemError> {

    // parse dependending on physics type
    match phys_transient {
        Some(TransientPhysics::Diff(phys)) => {
            parse_diff(parts, phys, dom0d_map, dom1d_map, scl0d_map, scl1d_map, var1d_map, line_num)?;
        },
        Some(TransientPhysics::DiffAdv(phys)) => {
            parse_diffadv(parts, phys, dom0d_map, dom1d_map, scl0d_map, scl1d_map, var1d_map, line_num)?;
        },
//        Some(TransientPhysics::DiffMulti(phys)) => {
//            parse_diff_multi(parts, phys, dom0d_map, dom1d_map, scl0d_map, scl1d_map, var1d_map, line_num)?;
//        },
//        Some(TransientPhysics::DiffAdvMulti(phys)) => {
//            parse_diffadv_multi(parts, phys, dom0d_map, dom1d_map, scl0d_map, scl1d_map, var1d_map, line_num)?;
//        },
        _ => {
            panic!("Line {}: Transient physics type not recognized.", line_num);
        },
    }

    // return
    Ok(())

}

fn parse_diff(
    parts: &Vec<String>,
    phys: &mut TransientDiff,
    dom0d_map: &std::collections::HashMap<String, usize>,
    dom1d_map: &std::collections::HashMap<String, usize>,
    scl0d_map: &std::collections::HashMap<String, usize>,
    scl1d_map: &std::collections::HashMap<String, usize>,
    var1d_map: &std::collections::HashMap<String, usize>,
    line_num: usize,
) -> Result<(), FVChemError> {

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
            let bnd_type = parts.get(1).expect(format!("Line {}: bnd requires 'bnd <bnd_type> <dom> ...'. Could not find bnd_type.", line_num).as_str());
            let dom_name = parts.get(2).expect(format!("Line {}: bnd requires 'bnd <bnd_type> <dom> ...'. Could not find dom.", line_num).as_str());
            let dom_id = *dom0d_map.get(dom_name).expect(format!("Line {}: dom0d name '{}' not found.", line_num, dom_name).as_str());

            // check boundary type
            match bnd_type.as_str() {
                "conc" => {
                    let c_name = parts.get(3).expect(format!("Line {}: conc requires 'bnd conc <dom> <conc>'. Could not find conc.", line_num).as_str());
                    let c_id = *scl0d_map.get(c_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, c_name).as_str());
                    phys.add_boundary_concentration(dom_id, c_id);
                },
                "flux" => {
                    let n_name = parts.get(3).expect(format!("Line {}: flux requires 'bnd flux <dom> <flux>'. Could not find flux.", line_num).as_str());
                    let n_id = *scl0d_map.get(n_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, n_name).as_str());
                    phys.add_boundary_flux(dom_id, n_id);
                },
                "mass_trn" => {
                    let k_name = parts.get(3).expect(format!("Line {}: mass_trn requires 'bnd mass_trn <dom> <masss_trn> <conc_ref>'. Could not find mass_trn.", line_num).as_str());
                    let c_name = parts.get(4).expect(format!("Line {}: mass_trn requires 'bnd mass_trn <dom> <masss_trn> <conc_ref>'. Could not find conc_ref.", line_num).as_str());
                    let k_id = *scl0d_map.get(k_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, k_name).as_str());
                    let c_id = *scl0d_map.get(c_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, c_name).as_str());
                    phys.add_boundary_masstransfer(dom_id, k_id, c_id);
                },
                _ => {
                    panic!("Line {}: Boundary type '{}' not recognized.", line_num, bnd_type);
                },
            }

        },
        "itr" => {

            // get interface type
            let itr_type = parts.get(1).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find itr_type.", line_num).as_str());
            let dom_a_name = parts.get(2).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find dom_a.", line_num).as_str());
            let dom_b_name = parts.get(3).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find dom_b.", line_num).as_str());
            let dom_a_id = *dom1d_map.get(dom_a_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_a_name).as_str());
            let dom_b_id = *dom1d_map.get(dom_b_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_b_name).as_str());

            // check interface type
            match itr_type.as_str() {
                "cont" => {
                    phys.add_interface_continuity(dom_a_id, dom_b_id);
                }
                "mass_trn" => {
                    let k_name = parts.get(4).expect(format!("Line {}: mass_trn requires 'itr mass_trn <dom_a> <dom_b> <masss_trn>'. Could not find mass_trn.", line_num).as_str());
                    let k_id = *scl1d_map.get(k_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, k_name).as_str());
                    phys.add_interface_masstransfer(dom_a_id, dom_b_id, k_id);
                }
                _ => {
                    panic!("Line {}: Interface type '{}' not recognized.", line_num, itr_type);
                }
            }

        },
        _ => {
            panic!("Line {}: Expected 'int', 'bnd', or 'itr' for internal, boundary, or interface condition.", line_num);
        },
    }

    // return
    Ok(())

}

fn parse_diffadv(
    parts: &Vec<String>,
    phys: &mut TransientDiffAdv,
    dom0d_map: &std::collections::HashMap<String, usize>,
    dom1d_map: &std::collections::HashMap<String, usize>,
    scl0d_map: &std::collections::HashMap<String, usize>,
    scl1d_map: &std::collections::HashMap<String, usize>,
    var1d_map: &std::collections::HashMap<String, usize>,
    line_num: usize,
) -> Result<(), FVChemError> {

    // get inputs
    let cond_type = parts.get(0).expect(format!("Line {}: Expected 'int', 'bnd', or 'itr' for internal, boundary, or interface condition.", line_num).as_str());

    // check if internal, boundary, or interface
    match cond_type.as_str() {
        "int" => {
            // get names of properties
            let dom_name = parts.get(1).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <velx> <src>'. Could not find dom.", line_num).as_str());
            let c_name = parts.get(2).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <velx> <src>'. Could not find conc.", line_num).as_str());
            let d_name = parts.get(3).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <velx> <src>'. Could not find diff.", line_num).as_str());
            let u_name = parts.get(4).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <velx> <src>'. Could not find velx.", line_num).as_str());
            let r_name = parts.get(5).expect(format!("Line {}: int requires 'int <dom> <conc> <diff> <velx> <src>'. Could not find src.", line_num).as_str());
            
            // get ids
            let dom_id = *dom1d_map.get(dom_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_name).as_str());
            let c_id = *var1d_map.get(c_name).expect(format!("Line {}: var1d name '{}' not found.", line_num, c_name).as_str());
            let d_id = *scl1d_map.get(d_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, d_name).as_str());
            let u_id = *scl1d_map.get(u_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, u_name).as_str());
            let r_id = *scl1d_map.get(r_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, r_name).as_str());
        
            // add internal domain to physics
            phys.add_domain(dom_id, c_id, d_id, u_id, r_id);
        },
        "bnd" => {
            // get boundary type
            let bnd_type = parts.get(1).expect(format!("Line {}: bnd requires 'bnd <bnd_type> <dom> ...'. Could not find bnd_type.", line_num).as_str());
            let dom_name = parts.get(2).expect(format!("Line {}: bnd requires 'bnd <bnd_type> <dom> ...'. Could not find dom.", line_num).as_str());
            let dom_id = *dom0d_map.get(dom_name).expect(format!("Line {}: dom0d name '{}' not found.", line_num, dom_name).as_str());

            // check boundary type
            match bnd_type.as_str() {
                "conc" => {
                    let c_name = parts.get(3).expect(format!("Line {}: conc requires 'bnd conc <dom> <conc>'. Could not find conc.", line_num).as_str());
                    let c_id = *scl0d_map.get(c_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, c_name).as_str());
                    phys.add_boundary_concentration(dom_id, c_id);
                },
                "flux" => {
                    let n_name = parts.get(3).expect(format!("Line {}: flux requires 'bnd flux <dom> <flux>'. Could not find flux.", line_num).as_str());
                    let n_id = *scl0d_map.get(n_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, n_name).as_str());
                    phys.add_boundary_flux(dom_id, n_id);
                },
                "flux_noadv" => {
                    let n_name = parts.get(3).expect(format!("Line {}: flux_noadv requires 'bnd flux_noadv <dom> <flux>'. Could not find flux.", line_num).as_str());
                    let n_id = *scl0d_map.get(n_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, n_name).as_str());
                    phys.add_boundary_fluxnoadvection(dom_id, n_id);
                }
                "mass_trn" => {
                    let k_name = parts.get(3).expect(format!("Line {}: mass_trn requires 'bnd mass_trn <dom> <masss_trn> <conc_ref>'. Could not find mass_trn.", line_num).as_str());
                    let c_name = parts.get(4).expect(format!("Line {}: mass_trn requires 'bnd mass_trn <dom> <masss_trn> <conc_ref>'. Could not find conc_ref.", line_num).as_str());
                    let k_id = *scl0d_map.get(k_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, k_name).as_str());
                    let c_id = *scl0d_map.get(c_name).expect(format!("Line {}: scl0d name '{}' not found.", line_num, c_name).as_str());
                    phys.add_boundary_masstransfer(dom_id, k_id, c_id);
                }
                _ => {
                    panic!("Line {}: Boundary type '{}' not recognized.", line_num, bnd_type);
                },
            }
        },
        "itr" => {
            // get interface type
            let itr_type = parts.get(1).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find itr_type.", line_num).as_str());
            let dom_a_name = parts.get(2).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find dom_a.", line_num).as_str());
            let dom_b_name = parts.get(3).expect(format!("Line {}: itr requires 'itr <itr_type> <dom_a> <dom_b> ...'. Could not find dom_b.", line_num).as_str());
            let dom_a_id = *dom1d_map.get(dom_a_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_a_name).as_str());
            let dom_b_id = *dom1d_map.get(dom_b_name).expect(format!("Line {}: dom1d name '{}' not found.", line_num, dom_b_name).as_str());

            // check interface type
            match itr_type.as_str() {
                "cont" => {
                    phys.add_interface_continuity(dom_a_id, dom_b_id);
                }
                "mass_trn" => {
                    let k_name = parts.get(4).expect(format!("Line {}: mass_trn requires 'itr mass_trn <dom_a> <dom_b> <masss_trn>'. Could not find mass_trn.", line_num).as_str());
                    let k_id = *scl1d_map.get(k_name).expect(format!("Line {}: scl1d name '{}' not found.", line_num, k_name).as_str());
                    phys.add_interface_masstransfer(dom_a_id, dom_b_id, k_id);
                }
                _ => {
                    panic!("Line {}: Interface type '{}' not recognized.", line_num, itr_type);
                }
            }
        },
        _ => {
            panic!("Line {}: Expected 'int', 'bnd', or 'itr' for internal, boundary, or interface condition.", line_num);
        },
    }

    // return
    Ok(())

}
