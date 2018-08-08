#![allow(non_snake_case)]
#![allow(unused_variables)]
extern crate ndarray;
extern crate rand;
use atoms::rand::prelude::*;
use atoms::ndarray::prelude::*;
use std::collections::HashMap;
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
    let mut ret: Vec<(usize, usize)> = Vec::new();
    let len = mol.dim().0;
    for i in 0..len {
        for j in 0.. len {
            if j <= i { continue; }
            for _ in 0..mol[[i, j]] {
                ret.push((i, j));
            }
        }
    }
    println!("edges");
    println!("{:?}", ret);
    ret
}
// need to bring chem shift in tables

// detects if rings are present in molecule
// returns None if no rings are present
// -returns vector containing carbons in ring if cycle is present
pub fn rings_present(mol: &Molecule) -> Option<Vec<(usize, usize)>> {
    let len = mol.dim().0;
    let mut nodes: Vec<usize> = Vec::new();
    for i in 0..len { // TODO find better way to initialize this
        nodes.push(i as usize);
    }
    let mut degree: Vec<usize> = vec![0; len];
    let mut singletons: Vec<usize> = Vec::new();
    let mut edge = edges(mol);
    while !nodes.is_empty() {
        //find degree of each nodes
        for (node1, node2) in edge.iter() {
            degree[*node1] += 1;
            degree[*node2] += 1;
        }
        // remove edges that contain nodes with degree of one
        if !degree.contains(&1) {
            // if there are no nodes with degree of 1 then the remaining nodes form a cycle
            // should return edges
            println!("Cycle detected");
            return Some(edge);
        }
        for od in degree.iter().enumerate() { // TODO - can probably make this for loop to an iter
            if *od.1 == 1 {
                singletons.push(od.0);
            }
        }
        nodes.retain(|node1| !singletons.contains(node1));
        edge.retain(|(node1,node2)| !singletons.contains(node1) && !singletons.contains(node2));
        degree = vec![0; len];
    }
    None
}
#[test]
fn test_rings_present() {
    let chrom1: Chromosome = vec![1,0,0,0,1,0,0,1,1,0];
    let test1 = chromosome_to_molecule(&chrom1, 5);
    let chrom2: Chromosome = vec![1,0,0,1,1,1,1,0,0,0,0,0,0,0,1];
    let test2 = chromosome_to_molecule(&chrom2, 6);
    assert_eq!(rings_present(&test1), None);
    assert_eq!(rings_present(&test2), Some(vec![(0,4),(0,5),(4,5)]));
}

pub fn compute_shift(
    mol: &Molecule, atoms: &Vec<&str>,
    chemical_formula: &HashMap<&str, i32>
) -> Vec<f32> {
    // find longest carbon chain
    let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
    let reduced_mol = reduce_molecule(mol);
    let r = rings_present(&reduced_mol);


    let temp: Vec<f32> = vec![0.0];
    temp
}
