#[macro_use(s)]
extern crate ndarray;
use  chem_shift::ndarray::prelude::*;
mod atoms;
use chemshift::atoms::*;

// need to bring  in tables

// detects if rings are present in molecule
// returns None if no rings are present
//returns vector containing carbons in ring if cycle is present
pub fn rings_present(mol: &Molecule, num_carbons: usize) => Option<Vec<Vec<usize>>> {
    let len = mol.dim().0;
    let mut subsets = Vec::new();
    let mut rings = Vec::new();
    // remove reflexive edges
    let mut edges = edges(mol).truncate(edges.len()/2);
    for (x, y) in edge {
        unimplemented();
        }
    }
    if !rings.is_empty() {
        Some(rings)
    } else {
        None
    }
}
pub fn compute_shift(
    mol: &Molecule, atoms: &Vec<&str>,
    chemical_formula: &HashMap<&str, i32>
) -> Vec<f32> {
    // find longest carbon chain
    let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
    let reduced_mol = reduce_molecule(mol);
    match rings_present(&reduced_mol, num_carbons) {
        Some(x) => (),
        None => (),
    }

}
