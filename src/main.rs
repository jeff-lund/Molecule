// Copyright (c) 2018 Jeff Lund
#![allow(unused_mut)]
use std::collections::HashMap;
use std::collections::HashSet;
use std::env::args;
use std::fs::File;
use std::io::prelude::*;
#[allow(dead_code)]
enum Atom {
    Hydrogen,
    DoubleBond,
    TripleBond,
    CH3,
    CH2,
    CH,
    CarboxylicAcid,
    Ester,
    Aldehyde,
    Ketone,
    Amide,
    Aromatic,
    Alkene,
    Alkyne,
    CO,
    CN,
    CCl,
    CBr,
}
#[allow(dead_code)]
struct Molecule {
    kind: Atom,
    bonds: Vec<i32>,
}
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
fn compute_ihd(elements: &HashMap<String, i32>) -> i32 {
    let c = match elements.get("C") {
        Some(n) => *n,
        None => 0
    };
    let h = match elements.get("H") {
        Some(n) => *n,
        None => 0
    };
    let n = match elements.get("N") {
        Some(n) => *n,
        None => 0
    };
    let cl = match elements.get("Cl") {
        Some(n) => *n,
        None => 0
    };
    let br = match elements.get("Br") {
        Some(n) => *n,
        None => 0
    };
    println!("C: {} H: {} N: {} Cl: {} Br: {}", c, h, n, cl, br);
    (2*c + 2 + n - h - cl - br) / 2
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

    Ok(())
}
