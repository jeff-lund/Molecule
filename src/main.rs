// Copyright (c) 2018 Jeff Lund
#![allow(dead_code)]
#![allow(unused_assignments)]
#![allow(unused_variables)]
mod utility;
use utility::*;
mod atoms;
use atoms::*;

use std::env::args;
use std::fs::File;
use std::io::prelude::*;

const DEBUG: bool = false;
const MAX_GENERATIONS: u32 = 200;
fn pprint (molecule: &Molecule, atoms: &Vec<&str>, peaks: &Vec<f32>) {
    for a in atoms {
        print!("|{} ", a);
    }
    println!("|");
    for row in molecule.structure.genrows() {
        println!("{}", row);
    }
    println!("");
    println!("Fitness: {}", molecule.fitness);
    print!("Carbon Assignments: {:?}", molecule.kind[0]);
    for entry in molecule.kind.iter().skip(1) {
        print!(", {:?}", entry);
    }
    println!("");
    print!("Chemical Shifts: {}", molecule.chemical_shifts[0]);
    for entry in molecule.chemical_shifts.iter().skip(1) {
        print!(", {}", entry);
    }
    println!("");
    print!("Experimental Chemical Shifts: {}", peaks[0]);
    for entry in peaks.iter().skip(1) {
        print!(", {}", entry);
    }
    println!("");
}

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
    if DEBUG == false {
        // Generate random molecule population
        let mut population: Vec<Molecule> = Vec::new();
        let mut pop = 0;
        while pop < POPULATION {
            let molecule = create_test_molecule(&atoms, bonds);
            if connected(&molecule.structure) && check_bonds(&molecule.structure, &atoms) {
                population.push(molecule);
                pop += 1;
            }
        }
        // START evolution
        let mut best_molecule: Molecule = Molecule::new(Structure::zeros((1,1)));
        for gen in 0..MAX_GENERATIONS {
            println!("Generation {}", gen);
            // calculate chemical shifts
            // calculate molecules fitness
            for molecule in population.iter_mut() {
                molecule.assign_carbons(&atoms);
                molecule.compute_shifts(&atoms);
                molecule.fitness(&peaks);
            }
            // check break condition
            best_molecule = best(&population);
            if best_molecule.fitness == 0.0 {
                break;
            }
            // Create new population from parent generation
            population = generate_children(population, &atoms, bonds)
            // END evolution
        }
        println!("Best fit found");
        pprint(&best_molecule, &atoms, &peaks);
    }
    if DEBUG == true {
        // DEBUG PRINTING
        println!("********************DEBUG*************************************");
        println!("{:?}", atoms);
        println!("IHD {}", ihd);
        println!("{:?}", chemical_formula);
        println!("{:?}", peaks);
        if symmetrical_carbons(&chemical_formula, &peaks) {
            println!("Symmetric carbons present");
        }
        println!("Bonds Assigned: {}", bonds);
        println!("********************END DEBUG*********************************");
        // END DEBUG PRINTING
        // START DEBUG CREATION
        let mut tot = 0;
        let mut fcon = 0;
        let mut fbond = 0;
        let mut failed_both = 0;
        let mut pop = 0;
        let mut population = Vec::new();
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

        for i in 0..10 {
            population[i].assign_carbons(&atoms);
            println!("{:?}", population[i].kind);
        }
        // END DEBUG CREATION
    }
    Ok(())
}
