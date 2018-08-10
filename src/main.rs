// Copyright (c) 2018 Jeff Lund
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

const DEBUG: bool = true;
const POPULATION: u32 = 64;

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    // parse file info
    let chemical_formula = parse_elemental_analysis(b.next().expect("invalid file"));
    let peaks = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&chemical_formula);
    let bonds = get_bonds(&chemical_formula);
    let atoms = get_atoms(&chemical_formula);

    // Generate random molecule population
    let mut population: Vec<Molecule> = Vec::new();
    if DEBUG == true {
        let mut tot = 0;
        let mut fcon = 0;
        let mut fbond = 0;
        let mut failed_both = 0;
        let mut pop = 0;
        while pop < POPULATION {
            let m = create_test_molecule(&atoms, bonds);
            let mut bond_flag = 0;
            let mut con_flag = 0;
            tot += 1;
            if !connected(&m.structure) {
                fcon += 1;
                con_flag = 1;
            }
            if !check_bonds(&m.structure, &atoms) {
                fbond += 1;
                bond_flag = 1;
            }
            if con_flag == 0 && bond_flag == 0 {
                    population.push(m);
                    pop += 1;
            }
            if con_flag == 1 && bond_flag == 1 {
                failed_both += 1;
            }
        }
        println!("Added {} molecules to population. Total molecules created: {}", pop, tot);
        println!("Failed bond check: {}", fbond);
        println!("Failed connected: {}", fcon);
        println!("Failed both: {}", failed_both);

        for p in population.iter_mut() {
            p.assign_carbons(&atoms);
            println!("{:?}", p.kind);
        }
    }
    // START evolution
    // calculate chemical shifts

    // calculate molecules fitness

    // Create new population from parent generation


    // randomly mutate new children

    // END evolution
    if DEBUG == true {
        // DEBUG PRINTING REMOVE LATER
        println!("********************DEBUG*************************************");
        println!("{:?}", atoms);
        println!("IHD {}", ihd);
        println!("{:?}", chemical_formula);
        println!("{:?}", peaks);
        if symmetrical_carbons(&chemical_formula, &peaks) {
            println!("Symmetric carbons present");
        }
        println!("Bonds - Total: {} | Assigned: {}", bonds.0, bonds.1);
        println!("********************END DEBUG*********************************");
        // END DEBUG PRINTING`
    }
    Ok(())
}
