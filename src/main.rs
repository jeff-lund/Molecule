// Copyright (c) 2018 Jeff Lund
#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_imports)]

#[macro_use(s)]
extern crate ndarray;
extern crate rand;
use rand::prelude::*;
use ndarray::prelude::*;

use std::collections::HashMap;
use std::collections::HashSet;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;

type Molecule = Array2<u8>;
type Chromosome = Vec<u8>;
// parsing elemental analysis
// breaks on single elements in formule i.e (CH4), needs error checking for valid elements
fn parse_elemental_analysis(formula: &str) -> HashMap<String, i32> {
    let mut chemical_formula = HashMap::new();
    let valid_elements: HashSet<_> = ["C", "H", "O", "N", "Cl", "Br"].iter().cloned().collect();
    let mut elements: Vec<&str> = formula.split(char::is_numeric).collect();
    elements.retain(|e| e != &"");
    for e in &elements {
        if !valid_elements.contains(e) {
            panic!("Invalid symbol in chemical formula")
        }
    }
    let mut quantities: Vec<&str> = formula.split(char::is_alphabetic).collect();
    quantities.retain(|e| e != &"");
    if elements.len() != quantities.len() {
        panic!("elements and quantities do not match. Invalid chemical formula");
    }
    for (elem, quant) in elements.iter().zip(quantities.iter()) {
        chemical_formula.insert(
            elem.to_string(),
            quant.parse::<i32>().expect("Not a number"),
        );
    }
    chemical_formula
}

fn parse_peaks(peaks: &str) -> Vec<i32> {
    let mut ret: Vec<i32> = Vec::new();
    let buf = peaks.split(',');
    for p in buf {
        ret.push(p.trim().parse::<i32>().expect("peaks has a non digit"));
    }
    ret.sort_by(|x, y| y.cmp(x));
    ret
}
/// Computes the index of hydrogen deficiency  to find level of unsaturation
// 1 degree of unsaturation = 1 ring or double bond in final structure
// IHD = (2C + 2 + N - H - X) / 2 where X is halogens
fn compute_ihd(elements: &HashMap<String, i32>) -> i32 {
    let mut ihd: i32 = 2;
    for (k, v) in elements.iter() {
        match k.as_str() {
            "C" => ihd += 2 * v,
            "H" => ihd -= v,
            "O" => continue,
            "N" => ihd += v,
            "Cl" => ihd -= v,
            "Br" => ihd -= v,
            _ => panic!("Unrecognized element in chemical formula"),
        }
    }
    ihd/2
}
#[test]
fn test_ihd() {
    let caffeine: HashMap<String, i32> = [("C".to_string(), 8), ("H".to_string(), 10),
        ("N".to_string(), 4), ("O".to_string(), 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&caffeine), 6);
    let acetic_acid: HashMap<String, i32> = [("C".to_string(), 2), ("H".to_string(), 4),
        ("O".to_string(), 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&acetic_acid), 1);
}

fn symmetrical_carbons(chemical_formula: &HashMap<String, i32>, chemical_peaks: &Vec<i32>) -> bool {
    let ncarbons = chemical_formula.get("C").expect("No carbons present in formula");
    let length = chemical_peaks.len() as i32;
    if *ncarbons == length {
        false
    } else if *ncarbons > length {
        true
    } else {
        panic!("Symmetrical carbons: cannot have more peaks than carbons");
    }
}

// Creates vec with individual atoms in C - O - N - Cl - Br order.
// Hydrogen atoms are not included in this vector
// Rework ugly copy paste code
fn get_atoms(chemical_formula: &HashMap<String, i32>) -> Vec<&str> {
    let mut v: Vec<&str> = Vec::new();
    match chemical_formula.get("C") {
        Some(n) => { for _ in 0..*n { v.push("C"); }},
        None => (),
    }
    match chemical_formula.get("O") {
        Some(n) => { for _ in 0..*n { v.push("O"); }},
        None => (),
    }
    match chemical_formula.get("N") {
        Some(n) => { for _ in 0..*n { v.push("N");}},
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(n) => {for _ in 0..*n { v.push("Cl"); }},
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(n) => { for _ in 0..*n { v.push("Br"); }},
        None => (),
    }
    v
}
// Returns (total bonds, assigned bonds)
// When building matrices only heavy atoms are assigned, hydrogen is ignored
// The total bonds are needed to check the final structure has exactly the number of
// open bonds to fill in with hydrogens
// Total bonds = (4 * carbon + 2 * oxygen + 3 * nitrogen + hydrogen + halogens) / 2
// assigned bonds = total bonds - hydrogen
fn get_bonds(chemical_formula: &HashMap<String, i32>) -> (i32, i32) {
    let mut total_bonds = 0;
    match chemical_formula.get("C") {
        Some(n) => total_bonds += n * 4,
        None => (),
    }
    match chemical_formula.get("O") {
        Some(n) => total_bonds += n * 2,
        None => (),
    }
    match chemical_formula.get("N") {
        Some(n) => total_bonds += n * 3,
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(n) => total_bonds += n,
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(n) => total_bonds += n,
        None => (),
    }
    let h: i32 = match chemical_formula.get("H") {
        Some(n) => *n,
        None => 0,
    };
    total_bonds = (total_bonds + h) / 2;
    let assigned_bonds = total_bonds - h;
    println!("Total: {} | Assigned: {}", total_bonds, assigned_bonds);
    (total_bonds, assigned_bonds)
}
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
fn molecule_to_chromosome(M: Molecule, length: usize) -> Chromosome {
    let mut chromosome = Chromosome::new();
    for i in 0..length {
        for j in 0..length {
            if j <= i { continue; }
            chromosome.push(M[[i, j]]);
        }
    }
    chromosome
}
// Randomly creates sample molecule
fn create_test_molecule(atoms: &Vec<&str>, bonds: (i32, i32)) -> Molecule {
    let mut rng = thread_rng();
    let l = atoms.len() as usize;
    let num_bonds = bonds.1;
    let chromosome_length: usize = (l * l - l)/2;
    let mut chromosome: Chromosome = vec![0; chromosome_length];

    let mut x: usize;
    for _ in 0..num_bonds {
        x = rng.gen_range(0, chromosome_length);
        chromosome[x] += 1;
    }
    println!("Chromesome: {:?}", chromosome);
    let ret = chromosome_to_molecule(&chromosome, l);
    println!("Adj Matrix:\n{:?}", ret);
    ret
}

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    let chemical_formula = parse_elemental_analysis(b.next().expect("invalid file"));
    let chemical_peaks: Vec<i32> = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&chemical_formula);
    let bonds = get_bonds(&chemical_formula);
    let atoms = get_atoms(&chemical_formula);
    let m = create_test_molecule(&atoms, bonds);
    println!("{:?}", molecule_to_chromosome(m, atoms.len()));

    println!("{:?}", atoms);
    println!("IHD {}", ihd);
    println!("{:?}", chemical_formula);
    println!("{:?}", chemical_peaks);

    if symmetrical_carbons(&chemical_formula, &chemical_peaks) {
        println!("Symmetric carbons present");
    }

    Ok(())
}
