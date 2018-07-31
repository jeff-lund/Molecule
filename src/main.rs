// Copyright (c) 2018 Jeff Lund
use std::env::args;
use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;

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
    bonds: [i32; 4],
}

fn main() -> std::io::Result<()> {
    // Read file from std::args to buffer
    let f = args().nth(1).expect("No file argument given");
    let mut file = File::open(f)?;
    let mut buffer =  String::new();
    file.read_to_string(&mut buffer)?;
    let mut b = buffer.lines(); // change this name later
    // parsing elemental analysis
    // breaks on single elements in formule i.e (CH4), needs error checking for valid elements
    let chemical_formula = b.next().expect("invalid file");
    let mut elements: Vec<&str> = chemical_formula.split(char::is_numeric).collect();
    elements.retain(|e| e != &"");
    let mut quantities: Vec<&str> = chemical_formula.split(char::is_alphabetic).collect();
    quantities.retain(|e| e != &"");
    let mut elemental_analysis = HashMap::new();
    if elements.len() != quantities.len() {
        panic!("elements and quantities do not match. Invalid chemical formula");
    }
    for (elem, quant) in elements.iter().zip(quantities.iter()) {
        elemental_analysis.insert(elem, quant.parse::<i32>().expect("Not a number"));
    }

    // parsing chemical peaks
    let mut chemical_peaks: Vec<i32> = Vec::new();
    let peaks = b.next().expect("invalid file").split(',');
    for p in peaks {
        chemical_peaks.push(p.trim().parse::<i32>().expect("peaks has a non digit"));
    }
    println!("{:?}", elemental_analysis);
    println!("{:?}", chemical_peaks);


    Ok(())
}
