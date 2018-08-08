#![allow(non_snake_case)]

#[derive(Debug, PartialEq)]
pub enum Atom {
    Hydrogen,
    CH3,            // 10-15
    CH2,            // 16-25
    CH,             // 20-50
    CCl,            // 10-65
    CBr,            // 10-65
    CN,             // 30-65
    CO,             // 50-90
    CNO,            // 50-90
    Alkene,         // 115-140
    Aromatic,       // 125-150
    Amide,          // 150-180
    CarboxylicAcid, // 160-185
    Ester,          // 160-185
    Aldehyde,       // 190-200
    Ketone,         // 205-220
}
#[derive(Debug)]
pub struct Molecule {
    kind: Atom,
    chemical_shift: i32,
}

impl Molecule {
    pub fn new(kind: Atom, shift: i32) -> Self {
        Molecule {kind, chemical_shift: shift}
    }
}
