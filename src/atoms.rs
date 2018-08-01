#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_mut)]

#[derive(Debug)]
pub enum Atom {
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
}

impl Molecule {
    pub fn new(kind: Atom, preset_bonds: i32) -> Self {
        let mut bonds = Vec::new();
        Molecule {kind, bonds, preset_bonds}
    }
    pub fn CH3() -> Molecule {
        Molecule::new(Atom::CH3, 3)
    }
    pub fn CH2() -> Molecule {
        Molecule::new(Atom::CH2, 2)
    }
    pub fn CH() -> Molecule {
        Molecule::new(Atom::CH3, 1)
    }
    pub fn CarboxylicAcid() -> Molecule {
        Molecule::new(Atom::CarboxylicAcid, 3)
    }
    pub fn Ester() -> Molecule {
        Molecule::new(Atom::Ester, 1)
    }
    pub fn Aldehyde() -> Molecule {
        Molecule::new(Atom::Aldehyde, 3)
    }
    pub fn Ketone() -> Molecule {
        Molecule::new(Atom::Ketone, 2)
    }
    pub fn Amide() -> Molecule {
        Molecule::new(Atom::Amide, 1)
    }
    pub fn Alkene() -> Molecule {
        Molecule::new(Atom::Alkene, 1)
    }
    pub fn Aromatic() -> Molecule {
        Molecule::new(Atom::Aromatic, 1)
    }
    pub fn CO() -> Molecule {
        Molecule::new(Atom::CO, 1)
    }
    pub fn COH() -> Molecule {
        Molecule::new(Atom::COH, 2)
    }
    pub fn CN() -> Molecule {
        Molecule::new(Atom::CN, 1)
    }
    pub fn CCl() -> Molecule {
        Molecule::new(Atom::CCl, 1)
    }
    pub fn CBr() -> Molecule {
        Molecule::new(Atom::CBr, 1)
    }
    //pub fn transform(&mut self, Atom) {}
}
