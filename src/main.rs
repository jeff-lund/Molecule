// Copyright (c) 2018 Jeff Lund
//#![allow(non_snake_case)]
#![allow(unused_imports)]
#![allow(dead_code)]
//#[macro_use(s)]
mod utility;
use utility::*;
mod atoms;
use atoms::*;
use std::collections::HashMap;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    // parse file info
    let chemical_formula = parse_elemental_analysis(b.next().expect("invalid file"));
    let chemical_peaks: Vec<f32> = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&chemical_formula);
    let bonds = get_bonds(&chemical_formula);
    let atoms = get_atoms(&chemical_formula);

    // Generate random molecule population
    let mut population: Vec<Molecule> = Vec::new();
    let mut pop = 0;
    while pop < 1000 {
        let m = create_test_molecule(&atoms, bonds);
        match m {
            Some(x) => {
                population.push(x);
                pop += 1;
            }
            None => (),
        }
    }
    println!("Added {} molecules to population", pop);
    // START evolution
    // calculate chemical shifts

    // calculate molecules fitness

    // Create new population from parent generation


    // randomly mutate new children

    // END evolution
    // DEBUG PRINTING REMOVE LATER
    println!("********************DEBUG*************************************");
    println!("{:?}", atoms);
    println!("IHD {}", ihd);
    println!("{:?}", chemical_formula);
    println!("{:?}", chemical_peaks);
    if symmetrical_carbons(&chemical_formula, &chemical_peaks) {
        println!("Symmetric carbons present");
    }
    println!("Bonds - Total: {} | Assigned: {}", bonds.0, bonds.1);
    println!("********************END DEBUG*********************************");
    // END DEBUG PRINTING

    Ok(())
}
