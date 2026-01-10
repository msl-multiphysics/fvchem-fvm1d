use std::fs::File;
use std::io::{self, BufRead, BufReader};

fn main() -> io::Result<()> {
    let file = File::open("examples/input_steady_file/r_face.csv")?;
    let reader = BufReader::new(file);

    let mut rows: Vec<(i32, f64, f64)> = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;

        // skip empty lines or header if needed
        if line.trim().is_empty() || i == 0 { continue; }

        let mut parts = line.split(',');

        let a: i32 = parts.next().unwrap().trim().parse().unwrap();
        let b: f64 = parts.next().unwrap().trim().parse().unwrap();
        let c: f64 = parts.next().unwrap().trim().parse().unwrap();

        rows.push((a, b, c));
    }

    println!("{rows:?}");
    Ok(())
}

// fn main() -> io::Result<()> {
//     let file = File::open("numbers.csv")?;
//     let reader = BufReader::new(file);
// 
//     let mut rows: Vec<(i32, i32, f64)> = Vec::new();
// 
//     for (i, line) in reader.lines().enumerate() {
//         let line = line?;
// 
//         // skip empty lines or header if needed
//         if line.trim().is_empty() || i == 0 { continue; }
// 
//         let mut parts = line.split(',');
// 
//         let a: i32 = parts.next().unwrap().trim().parse().unwrap();
//         let b: i32 = parts.next().unwrap().trim().parse().unwrap();
//         let c: f64 = parts.next().unwrap().trim().parse().unwrap();
// 
//         rows.push((a, b, c));
//     }
// 
//     println!("{rows:?}");
//     Ok(())
// }