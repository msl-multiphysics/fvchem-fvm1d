use fasteval::{Parser, Slab, Evaler};
use fvchem_fvm1d::*;
use shell_words;
use std::collections::HashMap;
use std::env;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader};
use std::sync::Arc;

fn main() -> Result<(), FVChemError> {
    
    // get path to solver and script
    let args: Vec<String> = env::args().collect();
    let path_scr = args[1].clone();

    // read file
    let file = File::open(&path_scr).expect("Failed to read script file.");
    let reader = BufReader::new(file);
    let lines = reader.lines();

    // initialize physics
    let mut phys_steady: Option<SteadyPhysics> = None;
    let mut phys_transient: Option<TransientPhysics> = None;
    let mut is_simu_defined = false;
    let mut is_phys_defined = false;
    let mut is_phys_steady = false;
    let mut phys_type: String;

    // initialize mesh
    let mut mesh = Mesh1D::default();
    let mut is_mesh_defined = false;

    // initialize problem
    let mut prob = Problem1D::default();
    let mut dom0d_map: HashMap<String, usize> = HashMap::new();
    let mut dom1d_map: HashMap<String, usize> = HashMap::new();
    let mut scl0d_map: HashMap<String, usize> = HashMap::new();
    let mut scl1d_map: HashMap<String, usize> = HashMap::new();
    let mut var1d_map: HashMap<String, usize> = HashMap::new();

    // initialize functions
    let mut func_map: HashMap<String, Arc<dyn Fn(usize, f64, &[f64]) -> f64 + Send + Sync>> = HashMap::new();
    let mut func_var: HashMap<String, Vec<usize>> = HashMap::new();

    // iterate through lines
    for (i, line) in lines.enumerate() {

        // parse line
        let line = line.expect(format!("Line {}: Failed to read line", i+1).as_str());
        let line = match line.split_once('#') {  // remove comments
            Some((s, _)) => s.trim().to_string(),
            None => line.trim().to_string(),
        };
        let parts = shell_words::split(&line).expect(format!("Line {}: Failed to parse line", i+1).as_str());

        // skip if empty
        if parts.len() == 0 {
            continue;
        }

        // run action depending on first word
        match parts[0].as_str() {
            "simu" => {
                
                // get inputs
                let time_type = parts.get(1).expect(format!("Line {}: simu requires 'simu <simu_type>'. Could not find simu_type.", i+1).as_str());

                // check inputs
                match time_type.as_str() {
                    "steady" => {
                        is_phys_steady = true;
                    },
                    "transient" => {
                        is_phys_steady = false;
                    },
                    _ => {
                        panic!("Line {}: Time_type must be 'steady' or 'transient'.", i+1);
                    }
                }

                // set physics defined flag
                is_simu_defined = true;

            },
            "mesh" => {

                // get inputs
                let mesh_type = parts.get(1).expect(format!("Line {}: mesh requires 'mesh <mesh_type> [...]'. Could not find mesh_type.", i+1).as_str());
                
                // create mesh depending on type
                match mesh_type.as_str() {
                    "bounds" => {  // mesh based on bounds
                        
                        // get inputs from args
                        let x_min: f64 = parts.get(2).expect(format!("Line {}: mesh bounds requires 'mesh bounds <x_min> <x_max> <num_cell>'. Could not find x_min.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse x_min into a number.", i+1).as_str());
                        let x_max: f64 = parts.get(3).expect(format!("Line {}: mesh bounds requires 'mesh bounds <x_min> <x_max> <num_cell>'. Could not find x_max.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse x_max into a number.", i+1).as_str());
                        let num_cell: usize = parts.get(4).expect(format!("Line {}: mesh bounds requires 'mesh bounds <x_min> <x_max> <num_cell>'. Could not find num_cell.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse num_cell into a number.", i+1).as_str());
                    
                        // create mesh
                        mesh = Mesh1D::new(x_min, x_max, num_cell)?;
                        is_mesh_defined = true;
                    
                    },
                    "file" => {  // mesh based on file
                        
                        // get inputs from args
                        let mesh_file = parts.get(2).expect(format!("Line {}: mesh file requires 'mesh file <mesh_path>'. Could not find mesh_path.", i+1).as_str());

                        // create mesh
                        mesh = Mesh1D::new_from_file(mesh_file.clone())?;
                        is_mesh_defined = true;

                    },
                    _ => {
                        panic!("Line {}: Unknown mesh_type '{}'.", i+1, mesh_type);
                    }
                }

            },
            "dom1d" => {
                
                // check that mesh is defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before domain", i+1)
                }

                // get inputs
                let name = parts.get(1).expect(format!("Line {}: dom1d requires 'dom1d <name> = <reg>'. Could not find name.", i+1).as_str());
                let eq_sign = parts.get(2).expect(format!("Line {}: dom1d requires 'dom1d <name> = <reg>'. Could not find '='.", i+1).as_str());
                let reg = parts.get(3).expect(format!("Line {}: dom1d requires 'dom1d <name> = <reg>'. Could not find reg.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse reg into a number.", i+1).as_str());

                // check inputs
                if name == "=" {
                    panic!("Line {}: dom1d name cannot be empty.", i+1);
                }
                if eq_sign != "=" {
                    panic!("Line {}: dom1d requires '=' after name.", i+1);
                }
                if dom1d_map.contains_key(name) {
                    panic!("Line {}: dom1d name '{}' already exists.", i+1, name);
                }

                // create domain
                let dom_id = prob.add_dom1d(&mesh, reg)?;
                dom1d_map.insert(name.to_string(), dom_id);

            },
            "dom0d" => {
                
                // check that mesh is defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before domain", i+1)
                }

                // get inputs
                let name = parts.get(1).expect(format!("Line {}: dom0d requires 'dom0d <name> = <dom1d_name> <bnd>'. Could not find name.", i+1).as_str());
                let eq_sign = parts.get(2).expect(format!("Line {}: dom0d requires 'dom0d <name> = <dom1d_name> <bnd>'. Could not find '='.", i+1).as_str());
                let dom1d_name = parts.get(3).expect(format!("Line {}: dom0d requires 'dom0d <name> = <dom1d_name> <bnd>'. Could not find dom1d_name.", i+1).as_str());
                let bnd = parts.get(4).expect(format!("Line {}: dom0d requires 'dom0d <name> = <dom1d_name> <bnd>'. Could not find bnd.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse bnd into a number.", i+1).as_str());

                // check inputs
                if name == "=" {
                    panic!("Line {}: dom0d name cannot be empty.", i+1);
                }
                if eq_sign != "=" {
                    panic!("Line {}: dom0d requires '=' after name.", i+1);
                }
                if dom0d_map.contains_key(name) {
                    panic!("Line {}: dom0d name '{}' already exists.", i+1, name);
                }

                // create domain
                let dom1d_id = match dom1d_map.get(dom1d_name) {
                    Some(id) => *id,
                    None => panic!("Line {}: dom1d name '{}' does not exist.", i+1, dom1d_name),
                };
                let dom_id = prob.add_dom0d(&mesh, dom1d_id, bnd)?;
                dom0d_map.insert(name.to_string(), dom_id);

            },
            "var1d" => {
                
                // check that mesh and simulation type are defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before variable", i+1)
                }
                if !is_simu_defined{
                    panic!("Line {}: simulation type must be defined before variable", i+1)
                }

                // get inputs 
                let name = parts.get(1).expect(format!("Line {}: var1d requires 'var1d <name> = <dom1d_name> <value_init> [...]'. Could not find name.", i+1).as_str());
                let eq_sign = parts.get(2).expect(format!("Line {}: var1d requires 'var1d <name> = <dom1d_name> <value_init> [...]'. Could not find '='.", i+1).as_str());
                let dom1d_name = parts.get(3).expect(format!("Line {}: var1d requires 'var1d <name> = <dom1d_name> <value_init> [...]'. Could not find dom1d_name.", i+1).as_str());
                let value_init = parts.get(4).expect(format!("Line {}: var1d requires 'var1d <name> = <dom1d_name> <value_init> [...]'. Could not find value_init.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse value_init into a number.", i+1).as_str());

                // check inputs
                if name == "=" {
                    panic!("Line {}: var1d name cannot be empty.", i+1);
                }
                if eq_sign != "=" {
                    panic!("Line {}: var1d requires '=' after name.", i+1);
                }
                if var1d_map.contains_key(name) {
                    panic!("Line {}: var1d name '{}' already exists.", i+1, name);
                }

                // create variable
                let dom1d_id = match dom1d_map.get(dom1d_name) {
                    Some(id) => *id,
                    None => panic!("Line {}: dom1d name '{}' does not exist.", i+1, dom1d_name),
                };
                let var_id = prob.add_var1d(dom1d_id, value_init)?;
                var1d_map.insert(name.to_string(), var_id);

                // get optional inputs
                let output_path = match parts.get(5) {
                    Some(path) => path.to_string(),
                    None => "".to_string(),
                };
                let output_step_str = match parts.get(6) {
                    Some(freq_str) => freq_str.to_string(),
                    None => "".to_string(),
                };

                // check optional inputs
                if is_phys_steady && output_path != "" {  // steady -> no output frequency
                    
                    // check output
                    if output_step_str != "" {
                        panic!("Line {}: var1d for steady physics with output requires 'var1d <name> = <dom1d_name> <value_init> <output_path>'. Found additional argument {}", i+1, output_step_str);
                    }

                    // set output
                    prob.set_var1d_write_steady(var_id, output_path.clone());

                    // create output directory
                    let output_vec = output_path.split('/').collect::<Vec<&str>>();
                    let output_folder = output_vec[..output_vec.len()-1].join("/");
                    create_dir_all(&output_folder).expect(format!("Line {}: Could not create output directory '{}'.", i+1, output_folder).as_str());

                } else if !is_phys_steady && output_path != "" {  // transient -> output frequency required
                    
                    // check output
                    if output_step_str == "" {
                        panic!("Line {}: var1d for transient physics with output requires 'var1d <name> = <dom1d_name> <value_init> <output_path> <output_step>'. Could not find output_step.", i+1);
                    }

                    // set output
                    let output_step: usize = output_step_str.parse().expect(format!("Line {}: Could not parse output_step into a number.", i+1).as_str());
                    prob.set_var1d_write_transient(var_id, output_path.clone(), output_step);

                    // create output directory
                    let output_vec = output_path.split('/').collect::<Vec<&str>>();
                    let output_folder = output_vec[..output_vec.len()-1].join("/");
                    create_dir_all(&output_folder).expect(format!("Line {}: Could not create output directory '{}'.", i+1, output_folder).as_str());

                }

            },
            "func" => {

                // get function name
                let name = parts.get(1).expect(format!("Line {}: func requires 'func <name> <var1_name> [...] = <expression>'. Could not find name.", i+1).as_str());

                // get variable names
                let eq_index = parts.iter().position(|x| x == "=").expect(format!("Line {}: func requires 'func <name> <var1_name> [...] = <expression>'. Could not find '='.", i+1).as_str());
                let var_name_vec = &parts[2..eq_index];
                
                // get variable ids and create variable index map
                let mut var_id_vec: Vec<usize> = Vec::new();
                let mut var_name_map: HashMap<String, usize> = HashMap::new();
                for (j, vn) in var_name_vec.iter().enumerate() { 
                    let var_id = match var1d_map.get(vn) {
                        Some(id) => *id,
                        None => panic!("Line {}: var1d name '{}' does not exist.", i+1, vn),
                    };
                    var_id_vec.push(var_id);
                    var_name_map.insert(vn.to_string(), j);
                }
                func_var.insert(name.to_string(), var_id_vec);

                // parse expression
                let expr_str = &parts[eq_index+1..].join(" ");
                let mut slab = Slab::new();
                let expr = Parser::new().parse(expr_str, &mut slab.ps).expect(format!("Line {}: Failed to parse expression.", i+1).as_str());

                // create function
                let func = Arc::new(move |t: usize, x: f64, vars: &[f64]| -> f64 {
                    let mut cb = |name: &str, _args: Vec<f64>| -> Option<f64> {
                        if name == "t" {
                            Some(t as f64)
                        } else if name == "x" {
                            Some(x)
                        } else if let Some(&i) = var_name_map.get(name) {
                            vars.get(i).copied()
                        } else {
                            panic!("Line {}: Unknown variable name '{}' in expression.", i+1, name);
                        }
                    };
                    expr.from(&slab.ps).eval(&slab, &mut cb).expect(format!("Line {}: Failed to evaluate expression.", i+1).as_str())
                });
                func_map.insert(name.to_string(), func);
                
            },
            "scl1d" => {

                // check that mesh and simulation type are defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before scalar", i+1)
                }
                if !is_simu_defined{
                    panic!("Line {}: simulation type must be defined before scalar", i+1)
                }

                // get inputs
                let name = parts.get(1).expect(format!("Line {}: scl1d requires 'scl1d <name> = <dom1d_name> <value> [...]'. Could not find name.", i+1).as_str());
                let eq_sign = parts.get(2).expect(format!("Line {}: scl1d requires 'scl1d <name> = <dom1d_name> <value> [...]'. Could not find '='.", i+1).as_str());
                let dom1d_name = parts.get(3).expect(format!("Line {}: scl1d requires 'scl1d <name> = <dom1d_name> <value> [...]'. Could not find dom1d_name.", i+1).as_str());
                let value_str = parts.get(4).expect(format!("Line {}: scl1d requires 'scl1d <name> = <dom1d_name> <value> [...]'. Could not find value.", i+1).as_str());

                // check inputs
                if name == "=" {
                    panic!("Line {}: scl1d name cannot be empty.", i+1);
                }
                if eq_sign != "=" {
                    panic!("Line {}: scl1d requires '=' after name.", i+1);
                }
                if scl1d_map.contains_key(name) {
                    panic!("Line {}: scl1d name '{}' already exists.", i+1, name);
                }

                // get domain 1d id
                let dom1d_id = match dom1d_map.get(dom1d_name) {
                    Some(id) => *id,
                    None => panic!("Line {}: dom1d name '{}' does not exist.", i+1, dom1d_name),
                };

                // create scalar depending on type
                match value_str.parse::<f64>().is_ok() {
                    true => {  // constant scalar
                        
                        // parse value
                        let value: f64 = value_str.parse().expect(format!("Line {}: Could not parse value into a number.", i+1).as_str());

                        // create scalar
                        let scl_id = prob.add_scl1d(dom1d_id, value)?;
                        scl1d_map.insert(name.to_string(), scl_id);

                    },
                    false => {
                        
                        // parse value
                        let func = match func_map.get(value_str) {
                            Some(f) => f,
                            None => panic!("Line {}: Function '{}' does not exist.", i+1, value_str),
                        };
                        let var_vec = match func_var.get(value_str) {
                            Some(v) => v.clone(),
                            None => panic!("Line {}: Function '{}' has no associated variables.", i+1, value_str),
                        };

                        // create scalar
                        let scl_id = prob.add_scl1d_from_function(dom1d_id, func.clone(), var_vec)?;
                        scl1d_map.insert(name.to_string(), scl_id);

                    }
                }

            },
            "scl0d" => {
                
                // check that mesh and simulation type are defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before scalar", i+1)
                }
                if !is_simu_defined{
                    panic!("Line {}: simulation type must be defined before scalar", i+1)
                }

                // get inputs
                let name = parts.get(1).expect(format!("Line {}: scl0d requires 'scl0d <name> = <dom0d_name> <value> [...]'. Could not find name.", i+1).as_str());
                let eq_sign = parts.get(2).expect(format!("Line {}: scl0d requires 'scl0d <name> = <dom0d_name> <value> [...]'. Could not find '='.", i+1).as_str());
                let dom0d_name = parts.get(3).expect(format!("Line {}: scl0d requires 'scl0d <name> = <dom0d_name> <value> [...]'. Could not find dom0d_name.", i+1).as_str());
                let value_str = parts.get(4).expect(format!("Line {}: scl0d requires 'scl0d <name> = <dom0d_name> <value> [...]'. Could not find value.", i+1).as_str());

                // check inputs
                if name == "=" {
                    panic!("Line {}: scl0d name cannot be empty.", i+1);
                }
                if eq_sign != "=" {
                    panic!("Line {}: scl0d requires '=' after name.", i+1);
                }
                if scl0d_map.contains_key(name) {
                    panic!("Line {}: scl0d name '{}' already exists.", i+1, name);
                }

                // get domain 0d id
                let dom0d_id = match dom0d_map.get(dom0d_name) {
                    Some(id) => *id,
                    None => panic!("Line {}: dom0d name '{}' does not exist.", i+1, dom0d_name),
                };

                // create scalar depending on type
                match value_str.parse::<f64>().is_ok() {
                    true => {  // constant scalar
                        // parse value
                        let value: f64 = value_str.parse().expect(format!("Line {}: Could not parse value into a number.", i+1).as_str());

                        // create scalar
                        let scl_id = prob.add_scl0d(dom0d_id, value)?;
                        scl0d_map.insert(name.to_string(), scl_id);
                    }
                    false => {

                        // parse value
                        let func = match func_map.get(value_str) {
                            Some(f) => f,
                            None => panic!("Line {}: Function '{}' does not exist.", i+1, value_str),
                        };
                        let var_vec = match func_var.get(value_str) {
                            Some(v) => v.clone(),
                            None => panic!("Line {}: Function '{}' has no associated variables.", i+1, value_str),
                        };

                        // create scalar
                        let scl_id = prob.add_scl0d_from_function(dom0d_id, func.clone(), var_vec)?;
                        scl0d_map.insert(name.to_string(), scl_id);

                    }
                }

            },
            "physics" => {
                
                // get inputs
                phys_type = parts.get(1).expect(format!("Line {}: physics requires 'physics <phys_type>'. Could not find phys_type.", i+1).as_str()).to_string();

                // set physics depending on type
                if is_phys_steady {
                    parse_steady_physics(&parts, &mut phys_steady, &phys_type, i+1)?;
                } else {
                    parse_transient_physics(&parts, &mut phys_transient, &phys_type, i+1)?;
                }
                
                // set physics defined flag
                is_phys_defined = true;

            },
            "int" | "bnd" | "itr" => {

                // check that mesh, simulation type, and physics are defined
                if !is_mesh_defined{
                    panic!("Line {}: mesh must be defined before {}", i+1, parts[0].as_str())
                }
                if !is_simu_defined{
                    panic!("Line {}: simulation type must be defined before {}", i+1, parts[0].as_str())
                }
                if !is_phys_defined{
                    panic!("Line {}: physics type must be defined before {}", i+1, parts[0].as_str())
                }

                // apply physics depend on time type
                if is_phys_steady {
                    parse_steady_condition(&parts, &mut phys_steady, &dom0d_map, &dom1d_map, &scl0d_map, &scl1d_map, &var1d_map, i+1)?;
                } else {
                    parse_transient_condition(&parts, &mut phys_transient, &dom0d_map, &dom1d_map, &scl0d_map, &scl1d_map, &var1d_map, i+1)?;
                }

            },
            "solve" => {

                // solve problem depending on time type
                if is_phys_steady {
                    // get inputs
                    let max_iter = parts.get(1).expect(format!("Line {}: solve for steady physics requires 'solve <max_iter> <tol> <damping>'. Could not find max_iter.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse max_iter into a number.", i+1).as_str());
                    let tol = parts.get(2).expect(format!("Line {}: solve for steady physics requires 'solve <max_iter> <tol> <damping>'. Could not find tolerance.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse tol into a number.", i+1).as_str());
                    let damping = parts.get(3).expect(format!("Line {}: solve for steady physics requires 'solve <max_iter> <tol> <damping>'. Could not find relaxation.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse damping into a number.", i+1).as_str());

                    // run solver
                    parse_steady_solver(&mut prob, &mut phys_steady, max_iter, tol, damping)?;
                } else {
                    // get inputs
                    let dt = parts.get(1).expect(format!("Line {}: solve for transient physics requires 'solve <dt> <num_step> <max_iter> <tol> <damping>'. Could not find dt.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse dt into a number.", i+1).as_str());
                    let num_step = parts.get(2).expect(format!("Line {}: solve for transient physics requires 'solve <dt> <num_step> <max_iter> <tol> <damping>'. Could not find num_step.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse num_step into a number.", i+1).as_str());
                    let max_iter = parts.get(3).expect(format!("Line {}: solve for transient physics requires 'solve <dt> <num_step> <max_iter> <tol> <damping>'. Could not find max_iter.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse max_iter into a number.", i+1).as_str());
                    let tol = parts.get(4).expect(format!("Line {}: solve for transient physics requires 'solve <dt> <num_step> <max_iter> <tol> <damping>'. Could not find tolerance.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse tol into a number.", i+1).as_str());
                    let damping = parts.get(5).expect(format!("Line {}: solve for transient physics requires 'solve <dt> <num_step> <max_iter> <tol> <damping>'. Could not find relaxation.", i+1).as_str()).parse().expect(format!("Line {}: Could not parse damping into a number.", i+1).as_str());
                    
                    // run solver
                    parse_transient_solver(&mut prob, &mut phys_transient, dt, num_step, max_iter, tol, damping)?;
                }

            },
            _ => {
                panic!("Unknown command '{}' in line {}", parts[0], i+1);
            }
        }

    }

    // return
    Ok(())

}
