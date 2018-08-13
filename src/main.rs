// Copyright (c) 2018 Jeff Lund
#![allow(unused_imports)]
#![allow(dead_code)]
mod utility;
use utility::*;
mod atoms;
use atoms::*;

use std::collections::HashMap;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;

const DEBUG: bool = true;
const DUMB_TESTING: bool = false;

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
    if symmetrical_carbons(&chemical_formula, &peaks) {
        eprintln!("There must be a chemical shift for each carbon atom in the chemical formula");
        std::process::exit(0);
    }
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

        //for p in population.iter_mut() {
        //    p.assign_carbons(&atoms);
    //        println!("{:?}", p.kind);
    //    }
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
    if DUMB_TESTING == true {
        // Figure out how to get this working as an integration test
        println!("test generate children");
        let t_atoms = vec!["C", "C","C","C","C","O","Cl"];
        let t_bonds = (16, 7);
        let mut t_population = Vec::new();
        let mut t_pop = 0;
        let mut t_runs = 0;
        let mut con_cnt = 0;
        let mut bnd_cnt = 0;
        while t_pop < POPULATION {
            let t_con: bool;
            let t_bnd: bool;
            t_runs += 1;
            if t_runs > 1000000 {
                println!("Failed bonds: {}", bnd_cnt);
                println!("Failed connection: {}", con_cnt);
                panic!("creation stuck in loop");
            }
            let mol = create_test_molecule(&t_atoms, t_bonds);
            t_bnd = check_bonds(&mol.structure, &t_atoms);
            t_con = connected(&mol.structure);
            if !t_bnd {
                bnd_cnt += 1;
            }
            if !t_con {
                con_cnt += 1;
            }
            if t_bnd && t_con {
                t_population.push(mol);
                t_pop += 1;
            }
        }
        let new_pop = generate_children(t_population, &t_atoms, t_bonds.1);
        assert_eq!(new_pop.len(), POPULATION);
        for p in new_pop {
            assert!(connected(&p.structure));
            assert!(check_bonds(&p.structure, &t_atoms));
        }
        println!("All good");
    }
    Ok(())
}
