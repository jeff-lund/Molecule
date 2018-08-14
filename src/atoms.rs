extern crate ndarray;
extern crate rand;
use atoms::rand::prelude::*;
use atoms::ndarray::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::cmp::Ordering::Less;
pub const POPULATION: usize = 64; // POPULATION must be even
pub const MUTATION_PROBABILITY: f64 = 0.05;
pub type Structure = Array2<u32>;
pub type Chromosome = Vec<u32>;

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
    // Aromatic,   -- Need better ring detection and maybe huckels rule fn to add this
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
    pub fitness: f32,
}


#[derive(Debug)]
pub struct Tree {
    pub alpha: Vec<usize>,
    pub beta: Vec<usize>,
    pub gamma: Vec<usize>,
    pub delta: Vec<usize>,
    pub epsilon: Vec<usize>,
}

// fitness init is ugly
impl Molecule {
    /// Creates new Molecule from a given structure. Kind and chemical shift vectors are initialized empty
    pub fn new(structure: Structure) -> Self {
        let kind: Vec<FunctionalGroup> = Vec::new();
        let chemical_shifts: Vec<f32> = Vec::new();
        Molecule {structure, kind, chemical_shifts, fitness: -999.9}
    }
    // Calculates fitness of a molecule by taking RMSD of chemical shifts of carbon atoms
    // fitness = sqrt( 1/N sum((chem_shift_calc - chem_shift_exp)^2))
    // Input: chemical peaks data
    pub fn fitness(&mut self, experimental: &Vec<f32>) {
        assert_eq!(self.chemical_shifts.len(), experimental.len(), "Peak data and chemical shifts should both represent the total number of carbon atoms in molecule.");
        let zipped = self.chemical_shifts.iter().zip(experimental.iter());
        self.fitness = (1.0/experimental.len() as f32) * (zipped.fold(0.0, |acc, (calc, exp)| acc + (calc-exp).powi(2))).sqrt();
    }
    /// Assigns functional groups to each carbon.
    // this needs to be broken up and cleaned up
    pub fn assign_carbons(&mut self, atoms: &Vec<&str>) {
        let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
        let edge_list = edges(&self.structure);
        for index in 0..num_carbons {
            let mut primary_edges = HashSet::new();
            let mut hydrogen_count;
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
            primary_edges.retain(|&(_x, y)| !dupes.contains(&y));
            hydrogen_count = 4 - primary_edges.len() - (dupes.len() * 2);
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
                     self.kind.push(Cyanide);
                     continue;
                 }
                if alkyne {
                    self.kind.push(Alkyne);
                    continue;
                }
                // match imine - could go more specific if needed
                if imine {
                    self.kind.push(Imine);
                    continue;
                }
                if alkene {
                    self.kind.push(Alkene);
                    continue;
                }
                // match carbonyl
                if carbonyl {
                    // carboxylic acid or ester or amide
                    if alcohol {
                        self.kind.push(CarboxylicAcid);
                        continue;
                    }
                    for (_a, b) in primary_edges.iter() {
                        if atoms[*b] == "O" {
                            self.kind.push(Ester);
                            continue;
                        } else if atoms[*b] == "N" {
                            self.kind.push(Amide);
                            continue;
                        }
                    }
                    // ketone or aldehyde or acyl halide
                    if hydrogen_count == 0 {
                        for (_a, b) in primary_edges.iter() {
                            if atoms[*b] == "Cl" {
                                self.kind.push(AcylChloride);
                                continue;
                            } else if atoms[*b] == "Br" {
                                self.kind.push(AcylBromide);
                                continue;
                            }
                        }
                        self.kind.push(Ketone);
                        continue;
                    } else {
                        self.kind.push(Aldehyde);
                        continue;
                    }
                }
            }
            // if only single bonds present
            // only bonded to carbons
            if secondary_atoms.iter().filter(|a| *a != &"C").count() == 0 {
                match hydrogen_count {
                    3 => self.kind.push(CH3),
                    2 => self.kind.push(CH2),
                    1 => self.kind.push(CH),
                    0 => self.kind.push(C),
                    _ => panic!("assign_bonds: Too many hydrogens!"),
                }
                continue;
            }
            if secondary_atoms.contains(&"Cl") {
                self.kind.push(CCl);
                continue;
            } else if secondary_atoms.contains(&"Br") {
                self.kind.push(CBr);
                continue;
            } else if secondary_atoms.contains(&"O") {
                self.kind.push(CO);
                continue;
            } else if secondary_atoms.contains(&"N") {
                self.kind.push(CN);
                continue;
            }
            panic!("No molecular assignment made");
        }
    }
    /// Computes chemical shift of each carbon in the structure
    /// See README for details about implemtation
    pub fn compute_shifts(&mut self, atoms: &Vec<&str>) {
        let num_carbons = atoms.iter().filter(|&c| *c == "C").count();
        for i  in 0..num_carbons {
            let mut tree = build_tree(i, &self.structure, atoms);

        }
    }
}
#[test]
fn test_assign_carbons() {
    let mut test = Molecule::new(chromosome_to_structure(&vec![1,2,1,0,0,0,1,0,0,0]));
    let atoms = vec!["C","C","O","O","N"];
    test.assign_carbons(&atoms);
    assert_eq!(test.kind, vec![CarboxylicAcid, CN]);
}

// Internal helper functions

/// Builds a tree of attached carbon atoms for each layer of chemical shift search
/// atoms refered to by index in atoms vec
// BUG removing edges will probably mess up cycles
fn build_tree(start: usize, structure: &Structure, atoms: &Vec<&str>)
-> Tree {
    let mut edges = edges(&structure);
    // START component vectors
    let alpha:       Vec<usize>;
    let mut beta:    Vec<usize> = Vec::new();
    let mut gamma:   Vec<usize> = Vec::new();
    let mut delta:   Vec<usize> = Vec::new();
    let mut epsilon: Vec<usize> = Vec::new();
    // END component vectors

    alpha = get_adjacent_nodes(start, &edges, atoms);
    // pull out any edge that contains the root
    edges.retain(|(x, y)| *x != start && *y != start);
    // beta layer
    for node in alpha.iter() {
        let mut temp = get_adjacent_nodes(*node, &edges, atoms);
        beta.append(&mut temp);
        edges.retain(|(x, y)| x != node && y != node);
    }
    // gamma layer
    for node in beta.iter() {
        let mut temp = get_adjacent_nodes(*node, &edges, atoms);
        gamma.append(&mut temp);
        edges.retain(|(x, y)| x != node && y != node);
    }
    // delta level
    for node in gamma.iter() {
        let mut temp = get_adjacent_nodes(*node, &edges, atoms);
        delta.append(&mut temp);
        edges.retain(|(x, y)| x != node && y != node);
    }
    //epsilon level
    for node in delta.iter() {
        let mut temp = get_adjacent_nodes(*node, &edges, atoms);
        epsilon.append(&mut temp);
        edges.retain(|(x, y)| x != node && y != node);
    }

    Tree { alpha, beta, gamma, delta, epsilon }
}
#[test]
fn test_build_tree() {
    let chrom = vec![1,1,0,0,0,0,0,0
                      ,0,1,1,0,0,0,0,
                         0,0,1,0,0,1,
                           0,0,0,0,0,
                             0,1,1,0,
                               0,0,0,
                                 0,0,
                                   0];
    let atoms = vec!["C","C","C","C","C","C","C","O","O"];
    println!("{:?}", chromosome_to_structure(&chrom));
    let test = build_tree(0, &chromosome_to_structure(&chrom), &atoms);
    let answer = Tree { alpha: vec![1, 2], beta: vec![4, 5],
        gamma: vec![6], delta: Vec::new(), epsilon: Vec::new() };
    assert_eq!(test.alpha, answer.alpha);
    assert_eq!(test.beta, answer.beta);
    assert_eq!(test.gamma, answer.gamma);
    assert_eq!(test.delta, answer.delta);
    assert_eq!(test.epsilon, answer.epsilon);
}
// gets all adjacent carbon nodes from a base node
fn get_adjacent_nodes(
    base_node: usize, edges: &Vec<(usize, usize)>, atoms: &Vec<&str>)
-> Vec<usize> {
    let mut ret = Vec::new();
    for (x, y) in edges {
        if *x == base_node && atoms[*y] == "C" && !ret.contains(y) {
            ret.push(*y);
        } else if *y == base_node && atoms[*x] == "C" && !ret.contains(x) {
            ret.push(*x);
        }
    }
    ret
}

/// Returns the length of a Structure.
/// As all Structures are square matrices either dimension can be returned
fn s_len(s: &Structure) -> usize {
    s.dim().0
}
/// Retruns the length of the molecule matrix derived from this chromosome
/// Derived from the quadratic equation
fn molecule_len(chromosome: &Chromosome) -> usize {
    ((1.0 + (1.0 + 8.0 * chromosome.len() as f32).sqrt()) / 2.0) as usize
}
#[test]
fn test_mol_len() {
    assert_eq!(molecule_len(&vec![0,0,0,0,0,0,0,0,0,0]), 5);
}
/// Returns (total bonds, assigned bonds)
/// When building matrices only heavy atoms are assigned, hydrogen is ignored
/// The total bonds are needed to check the final structure has exactly the number of
/// open bonds to fill in with hydrogens
/// Total bonds = (4 * carbon + 2 * oxygen + 3 * nitrogen + hydrogen + halogens) / 2
/// assigned bonds = total bonds - hydrogen
pub fn get_bonds(chemical_formula: &HashMap<&str, i32>) -> (u32, u32) {
    let mut total_bonds = 0;
    match chemical_formula.get("C") {
        Some(n) => total_bonds += *n as u32 * 4,
        None => (),
    }
    match chemical_formula.get("O") {
        Some(n) => total_bonds += *n as u32 * 2,
        None => (),
    }
    match chemical_formula.get("N") {
        Some(n) => total_bonds += *n as u32 * 3,
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(n) => total_bonds += *n as u32,
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(n) => total_bonds += *n as u32,
        None => (),
    }
    let h = match chemical_formula.get("H") {
        Some(n) => *n as u32,
        None => 0,
    };
    total_bonds = (total_bonds + h) / 2;
    let assigned_bonds = total_bonds - h;
    (total_bonds, assigned_bonds)
}
/// Transforms a chromsome into its matrix representation
pub fn chromosome_to_structure(chrom: &Chromosome) -> Structure {
    let length = molecule_len(chrom);
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
    let len = s_len(structure);
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
    let len = s_len(structure);
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
#[test]
fn test_connected() {
    // acyclic
    assert!(connected(&chromosome_to_structure(&vec![1,0,0,2,1,0])));
    // cyclic
    assert!(connected(&chromosome_to_structure(&vec![1,0,0,1,1,0,0,1,0,1])));
}
// checks that each atom in a molecule is not exceeding its maximum possible bonds
// C: 4 | O: 2 | N: 3 | Br: 1 | Cl: 1
// TODO try to enumerate rows and ditch row_index
pub fn check_bonds(structure: &Structure, atoms: &Vec<&str>) -> bool {
    let mut row_index = 0;
    let mut count;
    for row in structure.genrows() {
        count = row.iter().fold(0, |acc, x| acc + x);
        if count > match_element(atoms[row_index]) {
            return false;
        }
        row_index += 1;
    }
    true
}
#[test]
fn test_check_bonds_good() {
    let atoms = vec!["C", "O", "Cl"];
    let test = chromosome_to_structure(&vec![2,1,0]);
    assert!(check_bonds(&test, &atoms));

}
#[test]
#[should_panic]
fn test_check_bonds_panic() {
    let atoms = vec!["C", "O", "Cl"];
    let test = chromosome_to_structure(&vec![2,1,1]);
    assert!(check_bonds(&test, &atoms));
}

// Randomly creates sample molecule
// Needs tests for appropriate number of bonds and connectedness
// A chromosome is the concatenation of the values from the upper
// triangle exluding the diagonal
// TODO This algorithm is super inefficient find better way to generate graph
// most failures are in bonding checks
pub fn create_test_molecule(atoms: &Vec<&str>, bonds: (u32, u32)) -> Molecule {
    let mut rng = thread_rng();
    let l = atoms.len();
    let num_bonds = bonds.1;
    let chromosome_length: usize = (l * l - l)/2;
    let mut chromosome: Chromosome = vec![0; chromosome_length];

    let mut r: usize;
    for _ in 0..num_bonds {
        r = rng.gen_range(0, chromosome_length);
        // can't have more than a triple bond to any atom
        while chromosome[r] > 3 {
            r = rng.gen_range(0, chromosome_length);
        }
        chromosome[r] += 1;
    }
    Molecule::new(chromosome_to_structure(&chromosome))
}

/// Reduce moleculer structure to eliminate double/triple bonds
pub fn reduce_structure(structure: &Structure) -> Structure {
    let len = s_len(structure);
    let mut reduced = Structure::zeros((len, len));
    for i in 0..len {
        for j in 0..len {
            if structure[[i,j]] > 0 {
                reduced[[i, j]] = 1;
            }
        }
    }
    reduced
}
#[test]
fn test_reduce_structure() {
    assert_eq!(reduce_structure(&chromosome_to_structure(&vec![2,3,4,2,2,2,2,0,0,0])),
        chromosome_to_structure(&vec![1,1,1,1,1,1,1,0,0,0]));
}
/// Transforms an adjacency matrix into an edge list
/// The lower triangle is excluded as it is a mirror of the upper triangle
pub fn edges(structure: &Structure) -> (Vec<(usize, usize)>) {
    let mut ret: Vec<(usize, usize)> = Vec::new();
    let len = s_len(structure);
    for i in 0..len {
        for j in 0.. len {
            // only grab edges in the upper triangle
            if j <= i { continue; }
            for _ in 0..structure[[i, j]] {
                ret.push((i, j));
            }
        }
    }
    ret
}
#[test]
fn test_edges() {
    assert_eq!(edges(&chromosome_to_structure(&vec![1,1,1])), vec![(0,1), (0,2), (1,2)]);
}

/// Detects if rings are present in molecule
/// Returns None if no rings are present
/// Returns vector containing carbons in ring if cycle is present
/// or None if no rings(cycles) are present
/// Nodes contains nodes that are still possibilities for existing in a ring
/// Singletons contains nodes that are not in a ring
pub fn rings_present(structure: &Structure) -> Option<Vec<(usize, usize)>> {
    let len = s_len(structure);
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
    let test1 = chromosome_to_structure(&chrom1);
    let chrom2: Chromosome = vec![1,0,0,1,1,1,1,0,0,0,0,0,0,0,1];
    let test2 = chromosome_to_structure(&chrom2);
    assert_eq!(rings_present(&test1), None);
    assert_eq!(rings_present(&test2), Some(vec![(0,4),(0,5),(4,5)]));
}

/// Generates new child generation from parent population
/// Keeps best half of parent population and recombines parents to form children
/// Returns child population
pub fn generate_children(mut population:Vec<Molecule>, atoms: &Vec<&str>, num_bonds: u32) -> Vec<Molecule> {
    let mut children: Vec<Molecule> = Vec::new();
    let mut rng = thread_rng();
    // sort population
    population.sort_by(|a, b| (&a.fitness).partial_cmp(&b.fitness).unwrap_or(Less));
    // Take top half
    population.truncate(POPULATION/2);
    // Shuffle parents
    rng.shuffle(&mut population);
    // Create 2 children per pair
    while !population.is_empty() {
        let mol1 = population.pop().expect("odd number of molecules in population");
        let mol2 = population.pop().expect("odd number of molecules in population");
        let child1 = recombine(&mol1, &mol2, atoms, num_bonds);
        let child2 = recombine(&mol1, &mol2, atoms, num_bonds);
        children.push(child1);
        children.push(child2);
        children.push(mol1);
        children.push(mol2);
    }
    children
}
#[test]
fn test_generate_children() {
    let atoms = vec!["C", "C","C","C","C","O","Cl"];
    let bonds = (16, 7);
    let mut population = Vec::new();
    let mut pop = 0;
    let mut runs = 0;
    let mut con_cnt = 0;
    let mut bnd_cnt = 0;
    while t_pop < POPULATION {
        let t_con: bool;
        let t_bnd: bool;
        t_runs += 1;
        if t_runs > 1000000 {
            println!("Failed bonds: {}", bnd_cnt);
            println!("Failed connection: {}", con_cnt);
            panic!("creation stuck in loop");
        }
        let mol = create_test_molecule(&t_atoms, t_bonds);
        t_bnd = check_bonds(&mol.structure, &t_atoms);
        t_con = connected(&mol.structure);
        if !t_bnd {
            bnd_cnt += 1;
        }
        if !t_con {
            con_cnt += 1;
        }
        if t_bnd && t_con {
            t_population.push(mol);
            t_pop += 1;
        }
    }
    let new_pop = generate_children(t_population, &t_atoms, t_bonds.1);
    assert_eq!(new_pop.len(), POPULATION);
    for p in new_pop {
        assert!(connected(&p.structure));
        assert!(check_bonds(&p.structure, &t_atoms));
    }
}

/// recombine two parents to form a child chromosome
// TODO alter recombination probabilitoes based off of relative fitness
// TODO alter RP iteritively throughout generations
// TODO could send in RP and MP to iter alter thoughout gens
pub fn recombine(p1: &Molecule, p2: &Molecule, atoms: &Vec<&str>, num_bonds: u32) -> Molecule {
    let mut rng = thread_rng();
    let mut rand;
    let mut child: Chromosome = Vec::new();
    let chrom0 = structure_to_chromosome(&p1.structure);
    let chrom1 = structure_to_chromosome(&p2.structure);
    let mut runs = 0;
    loop {
        runs += 1;
        if runs > 1000000 {
            panic!("Recombine stuck in loop");
        }
        for i in 0..chrom1.len() {
            rand = rng.gen_range(0, 2);
            match rand {
                0 => child.push(chrom0[i]),
                1 => child.push(chrom1[i]),
                _ => panic!("recombine: generated rand out of range"),
            }
        }
        correct_child(&mut child, &chrom0, &chrom1, num_bonds);
        if rng.gen_bool(MUTATION_PROBABILITY) {
            mutate_child(&mut child, atoms, num_bonds)
        }
        if validate_chromosome(&child, atoms, num_bonds) {
            break;
        }
        child.clear();
    }
    Molecule::new(chromosome_to_structure(&child))
}
/// Mutates a random single bond in a molecule
pub fn mutate_child(chromosome: &mut Chromosome, atoms: &Vec<&str>, num_bonds: u32) {
    let mut rng = thread_rng();
    let l = chromosome.len();
    let mut gene_1;
    let mut gene_2;
    let mut runs= 0;
    loop {
        runs += 1;
        if runs > 10_000 {
            break;
        }
        gene_1 = rng.gen_range(0, l);
        gene_2 = rng.gen_range(0, l);
        if gene_1 != gene_2 {
            if chromosome[gene_1] == 0 { continue; }
            chromosome[gene_1] -= 1;
            chromosome[gene_2] += 1;
            if validate_chromosome(chromosome, atoms, num_bonds) {
                break;
            } else {
                chromosome[gene_1] += 1;
                chromosome[gene_2] -= 1;
            }
        }
    }
}
/// Validate the overall correctness of a chromosome
fn validate_chromosome(chromosome: &Chromosome, atoms: &Vec<&str>, num_bonds: u32) -> bool {
    // check proper number of bonds exist
    let bonds: u32 = chromosome.iter().sum();
    if bonds != num_bonds {
        return false;
    }
    // check to degree of each atom for bounds and connectedness of the moleculer structure
    let temp = chromosome_to_structure(chromosome);
    check_bonds(&temp, atoms) && connected(&temp)
}

/// Returns the maximum number of bonds a given element can have
pub fn match_element(element: &str) -> u32 {
    match element {
        "C" => 4,
        "O" => 2,
        "N" => 3,
        "Cl" | "Br" | "H" => 1,
        _ => panic!("unknown element found")
    }
}
/// Correct a chromosome such that it posesses the correct number of bonds
pub fn correct_child(child: &mut Chromosome, parent0: &Chromosome, parent1: &Chromosome, num_bonds: u32) {
    let mut rng = thread_rng();
    let mut r;
    let len = child.len();
    // correct number of bonds
    let mut current_bonds: u32 = child.iter().sum();
    let mut runs = 0;
    while current_bonds != num_bonds  {
        runs += 1;
        if runs > 100_000 {
            panic!("correct child stuck in loop");
        }
        // delete existing bond
        if current_bonds > num_bonds {
            r = rng.gen_range(0, len);
            while child[r] == 0 {
                r = rng.gen_range(0, len);
            }
            current_bonds -= child[r];
            child[r] = 0;
        }
        // add a bond from parent chromosome to child
        else if current_bonds < num_bonds {
            match rng.gen_range(0, 2) {
                0 => {
                    r = rng.gen_range(0, len);
                    while parent0[r] <= child[r] {
                        r = rng.gen_range(0, len);
                    }
                    current_bonds += parent0[r] - child[r];
                    child[r] = parent0[r];
                }
                1=> {
                    r = rng.gen_range(0, len);
                    while parent1[r] <= child[r] {
                        r = rng.gen_range(0, len);
                    }
                    current_bonds += parent1[r] - child[r];
                    child[r] = parent1[r];
                }
                _ => panic!("this isn't a thing"),
            }
        }
    }
}
