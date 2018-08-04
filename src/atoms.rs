#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_mut)]

#[derive(Debug)]
pub enum Atom {
    Hydrogen,
    CH3,            // 10-15
    CH2,            // 16-25
    CH,             // 20-50
    CCl,            // 10-65
    CBr,            // 10-65
    CN,             // 30-65
    CO,             // 50-90
    COH,            // 50-90
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
    bonds: Vec<i32>,
    preset_bonds: i32,
    chemical_shift: i32,
}

impl Molecule {
    pub fn new(kind: Atom, shift: i32, preset_bonds: i32) -> Self {
        let mut bonds = Vec::new();
        Molecule {kind, bonds, chemical_shift: shift, preset_bonds}
    }
    pub fn CH3(shift: i32) -> Molecule {
        Molecule::new(Atom::CH3, shift, 3)
    }
    pub fn CH2(shift: i32) -> Molecule {
        Molecule::new(Atom::CH2, shift, 2)
    }
    pub fn CH(shift: i32) -> Molecule {
        Molecule::new(Atom::CH3, shift, 1)
    }
    pub fn CarboxylicAcid(shift: i32) -> Molecule {
        Molecule::new(Atom::CarboxylicAcid, shift, 3)
    }
    pub fn Ester(shift: i32) -> Molecule {
        Molecule::new(Atom::Ester, shift, 1)
    }
    pub fn Aldehyde(shift: i32) -> Molecule {
        Molecule::new(Atom::Aldehyde, shift, 3)
    }
    pub fn Ketone(shift: i32) -> Molecule {
        Molecule::new(Atom::Ketone, shift, 2)
    }
    pub fn Amide(shift: i32) -> Molecule {
        Molecule::new(Atom::Amide, shift, 1)
    }
    pub fn Alkene(shift: i32) -> Molecule {
        Molecule::new(Atom::Alkene, shift, 1)
    }
    pub fn Aromatic(shift: i32) -> Molecule {
        Molecule::new(Atom::Aromatic, shift, 1)
    }
    pub fn CO(shift: i32) -> Molecule {
        Molecule::new(Atom::CO, shift, 1)
    }
    pub fn COH(shift: i32) -> Molecule {
        Molecule::new(Atom::COH, shift, 2)
    }
    pub fn CN(shift: i32) -> Molecule {
        Molecule::new(Atom::CN, shift, 1)
    }
    pub fn CCl(shift: i32) -> Molecule {
        Molecule::new(Atom::CCl, shift, 1)
    }
    pub fn CBr(shift: i32) -> Molecule {
        Molecule::new(Atom::CBr, shift, 1)
    }
    //pub fn transform(&mut self, Atom) {}
}
