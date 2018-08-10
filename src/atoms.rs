#![allow(non_snake_case)]
#![allow(unused_variables)]
extern crate ndarray;
extern crate rand;
use atoms::rand::prelude::*;
use atoms::ndarray::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;

pub type Structure = Array2<u8>;
pub type Chromosome = Vec<u8>;

#[derive(Debug, PartialEq)]
pub enum FunctionalGroup {
    CH3,
    CH2,
    CH,
    C,
    CCl,
    CBr,
    CN,
    CO,
    Alkene,
    Alkyne,
    Aromatic,
    Amide,
    Imine,
    CarboxylicAcid,
    Ester,
    Aldehyde,
    Ketone,
    AcylChloride,
    AcylBromide,
    Cyanide,
}
use FunctionalGroup::*;

#[derive(Debug)]
pub struct Molecule {
    pub structure: Structure,
    pub kind: Vec<FunctionalGroup>,
    pub chemical_shifts: Vec<f32>,
    pub fitness: Option<f32>,
}

impl Molecule {
    pub fn new(structure: Structure) -> Self {
        let kind: Vec<FunctionalGroup> = Vec::new();
        let chemical_shifts: Vec<f32> = Vec::new();
        Molecule {structure, kind, chemical_shifts, fitness: None}
    }
    // Calculates fitness of a molecule by taking RMSD of chemical shifts of carbon atoms
    // fitness = sqrt( 1/N sum((chem_shift_calc - chem_shift_exp)^2))
    // Input: chemical peaks data
    pub fn fitness(mut self, experimental: &Vec<f32>) {
        assert_eq!(self.chemical_shifts.len(), experimental.len(), "Peak data and chemical shifts should both represent the total number of carbon atoms in molecule.");
        let zipped = self.chemical_shifts.iter().zip(experimental.iter());
        self.fitness = Some((1.0/experimental.len() as f32) * (zipped.fold(0.0, |acc, (calc, exp)| acc + (calc-exp).powi(2))).sqrt());
    }
}

/// Assigns functional groups to each carbon.
// this needs to be broken up and cleaned up
pub fn assign_carbons(molecule: &mut Molecule, atoms: &Vec<&str>) {
    let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
    let edge_list = edges(&molecule.structure);
    for index in 0..num_carbons {
        let mut primary_edges = HashSet::new();
        let mut hydrogen_count: i32;
        let mut dupes = HashSet::new();
        let mut alcohol = false;
        // get primary edges
        for (n1, n2) in edge_list.iter() {
            if *n1 == index {
                if !primary_edges.insert((*n1, *n2)) {
                    dupes.insert(*n2);
                }
            } else if *n2 == index {
                if !primary_edges.insert((*n2, *n1)) {
                    dupes.insert(*n1);
                }
            }
        }
        // remove anything in dupes vec from primary edges, should only contain single bonded atoms
        primary_edges.retain(|&(x, y)| !dupes.contains(&y));
        hydrogen_count = 4 - primary_edges.len() as i32 - (dupes.len() * 2) as i32;
        if hydrogen_count < 0 {
            println!("dupes: {:?}", dupes);
            println!("primary_edges{:?}", primary_edges);
            println!("{:?}", molecule.structure);
            panic!("Negative hydrogen bonds");
        }
        let mut secondary_atoms = Vec::new();
        for (a, b) in primary_edges.iter() {
            secondary_atoms.push(atoms[*b]);
            if atoms[*b] == "O" {
                alcohol = edge_list.iter().filter(|(x, y)| (x != a && y == b) || (x == b && y != a)).count() == 0;
            }
        }
        // assign if double bonds present
        if !dupes.is_empty() {
            let mut carbonyl = false;
            let mut alkene = false;
            let mut alkyne = false;
            let mut imine = false;
            let mut amide = false;
            let mut cyanide = false;
            for d in dupes.iter() {
                match atoms[*d] {
                    "O" => carbonyl = true,
                    "N" => match imine {
                        true => cyanide = true,
                        false => imine = true,
                    }
                    "C" => match alkene {
                        true => alkyne = true,
                        false => alkene = true,
                    }
                    _ => panic!("unexpected element in duplicate matching"),
                }
            }
            // match triple bonds
            if cyanide {
                 molecule.kind.push(Cyanide);
                 continue;
             }
            if alkyne {
                molecule.kind.push(Alkyne);
                continue;
            }
            // match imine - could go more specific if needed
            if imine {
                molecule.kind.push(Imine);
                continue;
            }
            if alkene {
                molecule.kind.push(Alkene);
                continue;
            }
            // match carbonyl
            if carbonyl {
                // carboxylic acid or ester or amide
                if alcohol {
                    molecule.kind.push(CarboxylicAcid);
                    continue;
                }
                for (a, b) in primary_edges.iter() {
                    if atoms[*b] == "O" {
                        molecule.kind.push(Ester);
                        continue;
                    } else if atoms[*b] == "N" {
                        molecule.kind.push(Amide);
                        continue;
                    }
                }
                // ketone or aldehyde or acyl halide
                if hydrogen_count == 0 {
                    for (a, b) in primary_edges.iter() {
                        if atoms[*b] == "Cl" {
                            molecule.kind.push(AcylChloride);
                            continue;
                        } else if atoms[*b] == "Br" {
                            molecule.kind.push(AcylBromide);
                            continue;
                        }
                    }
                    molecule.kind.push(Ketone);
                    continue;
                } else {
                    molecule.kind.push(Aldehyde);
                    continue;
                }
            }
        }
        // if only single bonds present
        // only bonded to carbons
        if secondary_atoms.iter().filter(|a| *a != &"C").count() == 0 {
            match hydrogen_count {
                3 => molecule.kind.push(CH3),
                2 => molecule.kind.push(CH2),
                1 => molecule.kind.push(CH),
                0 => molecule.kind.push(C),
                _ => panic!("assign_bonds: Too many hydrogens!"),
            }
            continue;
        }
        if secondary_atoms.contains(&"Cl") {
            molecule.kind.push(CCl);
            continue;
        } else if secondary_atoms.contains(&"Br") {
            molecule.kind.push(CBr);
            continue;
        } else if secondary_atoms.contains(&"O") {
            molecule.kind.push(CO);
            continue;
        } else if secondary_atoms.contains(&"N") {
            molecule.kind.push(CN);
            continue;
        }
        panic!("No molecular assignment made")
    }
}

// length is the len x len dimensions of the desired matrix.
// Must be equal to the length of the atoms vector
pub fn chromosome_to_structure(chrom: &Chromosome, length: usize) -> Structure {
    let mut structure = Structure::zeros((length, length));
    let mut vec_index = 0;
    for i in 0..length {
        for j in 0..length {
            if  j <= i { continue; }
            structure[[i, j]] = chrom[vec_index];
            structure[[j, i]] = chrom[vec_index];
            vec_index += 1;
        }
    }
    structure
}
// creates a chromosome vector from a molecule matrix
pub fn structure_to_chromosome(structure: &Structure) -> Chromosome {
    let mut chromosome = Chromosome::new();
    let len = structure.dim().0;
    for i in 0..len {
        for j in 0..len {
            if j <= i { continue; }
            chromosome.push(structure[[i, j]]);
        }
    }
    chromosome
}
// Check if a moleculer structure is a connected graph using BFS
// components = 1 implies all vertices/atoms are in a single connected component
pub fn connected(structure: &Structure) -> bool {
    let mut components = 0;
    let len = structure.dim().0;
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
                    if marks[j] == 0 && structure[[v, j]] > 0 {
                        processing.push(j);
                    }
                }
            }
        }
    }
    components == 1
}

// checks that each atom in a molecule is not exceeding its maximum possible bonds
// C: 4 | O: 2 | N: 3 | Br: 1 | Cl: 1
// TODO try to enumerate rows and ditch row_index
pub fn check_bonds(structure: &Structure, atoms: &Vec<&str>) -> bool {
    let mut row_index = 0;
    let mut count;
    let mut allowed;
    for row in structure.genrows() {
        count = row.iter().fold(0, |acc, x| acc + x);
        allowed = match atoms[row_index] {
            "C" => 4,
            "N" => 3,
            "O" => 2,
            "H" | "Cl" | "Br" => 1,
            _ => panic!("check_bonds: unrecognized element in atoms array")
        };
        if count > allowed {
            return false;
        }
        row_index += 1;
    }
    true
}
#[test]
fn test_check_bonds_good() {
    let atoms = vec!["C", "O", "Cl"];
    let test = chromosome_to_structure(&vec![2,1,0], 3);
    assert!(check_bonds(&test, &atoms));

}
#[test]
#[should_panic]
fn test_check_bonds_panic() {
    let atoms = vec!["C", "O", "Cl"];
    let test = chromosome_to_structure(&vec![2,1,1], 3);
    assert!(check_bonds(&test, &atoms));
}

// Randomly creates sample molecule
// Needs tests for appropriate number of bonds and connectedness
// A chromosome is the concatenation of the values from the upper
// triangle exluding the diagonal
// TODO This algorithm is super inefficient find better way to generate graph
// most failures are in bonding checks
//pub fn create_test_molecule(atoms: &Vec<&str>, bonds: (i32, i32)) -> Option<Molecule> {
pub fn create_test_molecule(atoms: &Vec<&str>, bonds: (i32, i32)) -> Molecule {
    let mut rng = thread_rng();
    let l = atoms.len() as usize;
    let num_bonds = bonds.1;
    let chromosome_length: usize = (l * l - l)/2;
    let mut chromosome: Chromosome = vec![0; chromosome_length];

    let mut r: usize;
    for _ in 0..num_bonds {
        r = rng.gen_range(0, chromosome_length);
        // can't have more than a triple bond to any atom
        while chromosome[r] > 4 {
            r = rng.gen_range(0, chromosome_length);
        }
        chromosome[r] += 1;
    }
    let ret = Molecule::new(chromosome_to_structure(&chromosome, l));
    /*
    if check_bonds(&ret.structure, &atoms) && connected(&ret.structure) {
        Some(ret)
    } else {
        None
    }*/
    ret // for testing only
}

// reduce molecule to eliminate double/triple bonds
pub fn reduce_structure(structure: &Structure) -> Structure {
    let mut chrom = structure_to_chromosome(structure);
    for i in 0..chrom.len() {
        match chrom[i] {
            0 => (),
            _ => chrom[i] = 1,
        }
    }
    chromosome_to_structure(&chrom, structure.dim().0)
}

pub fn edges(structure: &Structure) -> (Vec<(usize, usize)>) {
    let mut ret: Vec<(usize, usize)> = Vec::new();
    let len = structure.dim().0;
    for i in 0..len {
        for j in 0.. len {
            if j <= i { continue; }
            for _ in 0..structure[[i, j]] {
                ret.push((i, j));
            }
        }
    }
    ret
}
// need to bring chem shift in tables

// detects if rings are present in molecule
// returns None if no rings are present
// -returns vector containing carbons in ring if cycle is present
// nodes contains nodes that are still possibilities for existing in a rings
// singletons contains nodes that are not in a ring
pub fn rings_present(structure: &Structure) -> Option<Vec<(usize, usize)>> {
    let len = structure.dim().0;
    let mut nodes: Vec<usize> = Vec::new();
    for i in 0..len { // TODO find better way to initialize this
        nodes.push(i as usize);
    }
    let mut degree: Vec<usize> = vec![0; len];
    let mut singletons: Vec<usize> = Vec::new();
    let mut edge_list = edges(&structure);

    while !nodes.is_empty() {
        //find degree of each nodes
        for (node1, node2) in edge_list.iter() {
            degree[*node1] += 1;
            degree[*node2] += 1;
        }
        if !degree.contains(&1) {
            return Some(edge_list);
        }
        // Add nodes with degree of 1 to singtons vec, theyre not in the ring
        for d in degree.iter().enumerate() { // TODO - can probably make this for loop to an iter
            if *d.1 == 1 {
                singletons.push(d.0);
            }
        }
        nodes.retain(|node1| !singletons.contains(node1));
        edge_list.retain(|(node1,node2)| !singletons.contains(node1) && !singletons.contains(node2));
        degree = vec![0; len];
    }
    None
}
#[test]
fn test_rings_present() {
    let chrom1: Chromosome = vec![1,0,0,0,1,0,0,1,1,0];
    let test1 = chromosome_to_structure(&chrom1, 5);
    let chrom2: Chromosome = vec![1,0,0,1,1,1,1,0,0,0,0,0,0,0,1];
    let test2 = chromosome_to_structure(&chrom2, 6);
    assert_eq!(rings_present(&test1), None);
    assert_eq!(rings_present(&test2), Some(vec![(0,4),(0,5),(4,5)]));
}

pub fn compute_shift(
    molecule: &Molecule, atoms: &Vec<&str>,
    chemical_formula: &HashMap<&str, i32>
) -> Vec<f32> {
    /*
    // find longest carbon chain
    let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
    let reduced_mol = reduce_molecule(molecule.structure);
    let r = rings_present(&reduced_mol);


    let temp: Vec<f32> = vec![0.0];
    temp */
    unimplemented!();
}
