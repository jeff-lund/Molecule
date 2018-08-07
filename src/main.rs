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
    let mut i: i32;
    match chemical_formula.get("C") {
        Some(i) => {
            for _ in 0..*i {
                v.push("C");
            }
        },
        None => (),
    }
    match chemical_formula.get("O") {
        Some(i) => {
            for _ in 0..*i {
                v.push("O");
            }
        },
        None => (),
    }
    match chemical_formula.get("N") {
        Some(i) => {
            for _ in 0..*i {
                v.push("N");
            }
        },
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(i) => {
            for _ in 0..*i {
                v.push("Cl");
            }
        },
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(i) => {
            for _ in 0..*i {
                v.push("Br");
            }
        },
        None => (),
    }
    v
}
fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    let mut chemical_formula = parse_elemental_analysis(b.next().expect("invalid file"));
    let mut chemical_peaks: Vec<i32> = parse_peaks(b.next().expect("invalid file"));
    let ihd = compute_ihd(&chemical_formula);
    let atoms = get_atoms(&chemical_formula);
    println!("{:?}", atoms);
    println!("IHD {}", ihd);
    println!("{:?}", chemical_formula);
    println!("{:?}", chemical_peaks);

    if symmetrical_carbons(&chemical_formula, &chemical_peaks) {
        println!("Symmetric carbons present");
    } else {
        println!("All carbons accounted for in peaks")
    }

    Ok(())
}
