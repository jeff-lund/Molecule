// Copyright (c) 2018 Jeff Lund
#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_mut)]
#![allow(unused_variables)]
#![allow(unused_imports)]
use std::collections::HashMap;
use std::collections::HashSet;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;
mod atoms;
use atoms::Atom;
use atoms::Molecule;

// parsing elemental analysis
// breaks on single elements in formule i.e (CH4), needs error checking for valid elements
fn parse_chemical_formula(formula: &str) -> HashMap<String, i32> {
    let mut elemental_analysis = HashMap::new();
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
        elemental_analysis.insert(
            elem.to_string(),
            quant.parse::<i32>().expect("Not a number"),
        );
    }
    elemental_analysis
}

fn parse_peaks(peaks: &str) -> Vec<i32> {
    let mut ret: Vec<i32> = Vec::new();
    let buf = peaks.split(',');
    for p in buf {
        ret.push(p.trim().parse::<i32>().expect("peaks has a non digit"));
    }
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
//add tests with halogens
fn test_ihd() {
    let caffeine: HashMap<String, i32> = [("C".to_string(), 8), ("H".to_string(), 10),
        ("N".to_string(), 4), ("O".to_string(), 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&caffeine), 6);
    let acetic_acid: HashMap<String, i32> = [("C".to_string(), 2), ("H".to_string(), 4),
        ("O".to_string(), 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&acetic_acid), 1);
}

fn build_elements(elemental_analysis: &HashMap<String, i32>, chemical_peaks: &Vec<i32>, ihd: i32) -> Vec<Molecule> {
    let mut build: Vec<Molecule> = Vec::new();

    let mut N_avail: i32 = match elemental_analysis.get("N") {
        Some(n) => *n as i32,
        None => 0,
    };
    let mut O_avail: i32 = match elemental_analysis.get("O") {
        Some(n) => *n as i32,
        None => 0,
    };
    for shift in chemical_peaks.iter() {
        if *shift <= 15 {
            build.push(Molecule::CH3());
        } else if *shift <= 25 {
            build.push(Molecule::CH2());
        } else if *shift <= 50 {
            build.push(Molecule::CH());
        } else if *shift <= 90 {
            if O_avail > 0 {
                build.push(Molecule::COH());
                O_avail -= 1;
            } else if N_avail > 0 {
                build.push(Molecule::CN());
                N_avail -= 1;
            } else { panic!("shift in the 50-90 with no O or N available")}
        } else if *shift <= 125 {
            build.push(Molecule::Alkene());
        } else if *shift <= 150 {
            build.push(Molecule::Aromatic());
        } else if *shift <= 170 {
            build.push(Molecule::Ester());
            O_avail -= 2;
        } else if *shift < 190 {
            build.push(Molecule::CarboxylicAcid());
            O_avail -= 2;
        } else if *shift <= 205 {
            build.push(Molecule::Aldehyde());
            O_avail -= 1;
        } else if *shift <= 220 {
            build.push(Molecule::Ketone());
            O_avail -= 1;
        } else {
            panic!("Chemical Shift out of range. Values must be  less than 220 cm-1");
        }
    }
    println!("{:?}", build);
    build
}

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    let mut elemental_analysis = parse_chemical_formula(b.next().expect("invalid file"));
    let mut chemical_peaks: Vec<i32> = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&elemental_analysis);
    println!("IHD {}", ihd);
    println!("{:?}", elemental_analysis);
    println!("{:?}", chemical_peaks);
    let mut molcules: Vec<Molecule> = build_elements(&elemental_analysis, &chemical_peaks, ihd);


    Ok(())
}
