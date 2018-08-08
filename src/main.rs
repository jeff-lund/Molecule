// Copyright (c) 2018 Jeff Lund
//#![allow(non_snake_case)]
#![allow(unused_imports)]
#![allow(dead_code)]
//#[macro_use(s)]
extern crate ndarray;
extern crate rand;
use rand::prelude::*;
use ndarray::prelude::*;

use std::collections::HashMap;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;

mod utility;
use utility::*;
type Molecule = Array2<u8>;
type Chromosome = Vec<u8>;

// length is the len x len dimensions of the desired matrix.
// Must be equal to the length of the atoms vector
fn chromosome_to_molecule(chrom: &Chromosome, length: usize) -> Molecule {
    let mut m = Molecule::zeros((length, length));
    let mut c = 0;
    for i in 0..length {
        for j in 0..length {
            if  j <= i { continue; }
            m[[i, j]] = chrom[c];
            m[[j, i]] = chrom[c];
            c += 1;
        }
    }
    m
}
// creates a chromosome vector from a molecule matrix
// length is the length of the atoms vector, aka len for len x len matrix
fn molecule_to_chromosome(mol: Molecule, length: usize) -> Chromosome {
    let mut chromosome = Chromosome::new();
    for i in 0..length {
        for j in 0..length {
            if j <= i { continue; }
            chromosome.push(mol[[i, j]]);
        }
    }
    chromosome
}
// Check if a Molecule is a connected graph using BFS
// components = 1 implies all vertices/atoms are in a single connected component
fn connected(molecule: &Molecule, len: usize) -> bool {
    let mut components = 0;
    let mut marks: Vec<usize> = vec![0; len];
    let mut processing: Vec<usize> = Vec::new();
    let mut v;
    for i in 0..len {
        if marks[i] == 0 {
            components += 1;
            processing.push(marks[i]);
            while !processing.is_empty() {
                v = processing.remove(0);
                marks[v] += components;
                for j in 0..len {
                    if marks[j] == 0 && molecule[[v, j]] > 0 {
                        processing.push(j);
                    }
                }
            }
        }
    }
    components == 1
}
// Randomly creates sample molecule
// Needs tests for appropriate number of bonds and connectedness
fn create_test_molecule(atoms: &Vec<&str>, bonds: (i32, i32)) -> Option<Molecule> {
    let mut rng = thread_rng();
    let l = atoms.len() as usize;
    let num_bonds = bonds.1;
    let chromosome_length: usize = (l * l - l)/2;
    let mut chromosome: Chromosome = vec![0; chromosome_length];

    let mut x: usize;
    for _ in 0..num_bonds {
        x = rng.gen_range(0, chromosome_length);
        while chromosome[x] > 4 {
            x = rng.gen_range(0, chromosome_length);
        }
        chromosome[x] += 1;
    }
    let ret = chromosome_to_molecule(&chromosome, l);
    if !connected(&ret, l) {
        None
    } else {
        Some(ret)
    }
}

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    let chemical_formula = parse_elemental_analysis(b.next().expect("invalid file"));
    let chemical_peaks: Vec<f32> = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&chemical_formula);
    let bonds = get_bonds(&chemical_formula);
    let atoms = get_atoms(&chemical_formula);

    // Generate random molecule population
    let mut population: Vec<Molecule> = Vec::new();
    let mut pop = 0;
    while pop < 100 {
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
