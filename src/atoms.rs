#![allow(non_snake_case)]
extern crate ndarray;
extern crate rand;
use atoms::rand::prelude::*;
use atoms::ndarray::prelude::*;

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
pub struct FuncGroup {
    kind: Atom,
    chemical_shift: i32,
}

impl FuncGroup {
    pub fn new(kind: Atom, shift: i32) -> Self {
        FuncGroup {kind, chemical_shift: shift}
    }
}

pub type Molecule = Array2<u8>;
pub type Chromosome = Vec<u8>;

// length is the len x len dimensions of the desired matrix.
// Must be equal to the length of the atoms vector
pub fn chromosome_to_molecule(chrom: &Chromosome, length: usize) -> Molecule {
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
pub fn molecule_to_chromosome(mol: &Molecule) -> Chromosome {
    let mut chromosome = Chromosome::new();
    let len = mol.dim().0;
    for i in 0..len {
        for j in 0..len {
            if j <= i { continue; }
            chromosome.push(mol[[i, j]]);
        }
    }
    chromosome
}
// Check if a Molecule is a connected graph using BFS
// components = 1 implies all vertices/atoms are in a single connected component
pub fn connected(molecule: &Molecule, len: usize) -> bool {
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
pub fn create_test_molecule(atoms: &Vec<&str>, bonds: (i32, i32)) -> Option<Molecule> {
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
    if connected(&ret, l) {
        Some(ret)
    } else {
        None
    }
}

// reduce molecule to eliminate double/truple bonds
pub fn reduce_molecule(mol: &Molecule) -> Molecule {
    let mut chrom = molecule_to_chromosome(mol);
    for i in 0..chrom.len() {
        match chrom[i] {
            0 => (),
            _ => chrom[i] = 1,
        }
    }
    chromosome_to_molecule(&chrom, mol.dim().0)
}

pub fn edges(mol: &Molecule) -> (Vec<(usize, usize)>) {
    let mut ret: Vec<(usize, usize)> = Vec::new()
    let len = mol.dim().0;
    for i in 0..len {
        for j in 0.. len {
            for 0..mol[[i, j]] {
                ret.push((i, j));
            }
        }
    }
    ret
}
