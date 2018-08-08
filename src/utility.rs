use std::collections::HashMap;
use std::collections::HashSet;

// Parses elemental analysis into a HashMap for easy access
// to the elements of the chemical formula
// Requires explicit counting for all atoms in formula, i.e. C2H6O1 not C2H6O
pub fn parse_elemental_analysis(formula: &str) -> HashMap<&str, i32> {
    let mut chemical_formula = HashMap::new();
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
        chemical_formula.insert(
            *elem,
            quant.parse::<i32>().expect("Not a number"),
        );
    }
    chemical_formula
}
// Collects input peak values into vector
pub fn parse_peaks(peaks: &str) -> Vec<f32> {
    let mut ret: Vec<f32> = Vec::new();
    let buf = peaks.split(',');
    for p in buf {
        ret.push(p.trim().parse::<f32>().expect("peaks incorrectly formatted, expects floats"));
    }
    ret.sort_by(|x, y| y.partial_cmp(x).unwrap());
    ret
}

pub fn symmetrical_carbons(chemical_formula: &HashMap<&str, i32>,
        chemical_peaks: &Vec<f32>) -> bool {
    let ncarbons = chemical_formula.get("C").expect("No carbons present in formula");
    let length = chemical_peaks.len() as i32;
    if *ncarbons == length {
        false
    } else if *ncarbons > length {
        true
    } else {
        panic!("Symmetrical carbons: cannot have more peaks than carbons");
    }
}

// Creates vec with individual atoms in C - O - N - Cl - Br order.
// Hydrogen atoms are not included in this vector
// Rework ugly copy paste code
pub fn get_atoms(chemical_formula: &HashMap<&str, i32>) -> Vec<&'static str> {
    let mut v: Vec<&str> = Vec::new();
    match chemical_formula.get("C") {
        Some(n) => { for _ in 0..*n { v.push("C"); }},
        None => (),
    }
    match chemical_formula.get("O") {
        Some(n) => { for _ in 0..*n { v.push("O"); }},
        None => (),
    }
    match chemical_formula.get("N") {
        Some(n) => { for _ in 0..*n { v.push("N");}},
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(n) => {for _ in 0..*n { v.push("Cl"); }},
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(n) => { for _ in 0..*n { v.push("Br"); }},
        None => (),
    }
    v
}
// Returns (total bonds, assigned bonds)
// When building matrices only heavy atoms are assigned, hydrogen is ignored
// The total bonds are needed to check the final structure has exactly the number of
// open bonds to fill in with hydrogens
// Total bonds = (4 * carbon + 2 * oxygen + 3 * nitrogen + hydrogen + halogens) / 2
// assigned bonds = total bonds - hydrogen
pub fn get_bonds(chemical_formula: &HashMap<&str, i32>) -> (i32, i32) {
    let mut total_bonds = 0;
    match chemical_formula.get("C") {
        Some(n) => total_bonds += n * 4,
        None => (),
    }
    match chemical_formula.get("O") {
        Some(n) => total_bonds += n * 2,
        None => (),
    }
    match chemical_formula.get("N") {
        Some(n) => total_bonds += n * 3,
        None => (),
    }
    match chemical_formula.get("Cl") {
        Some(n) => total_bonds += n,
        None => (),
    }
    match chemical_formula.get("Br") {
        Some(n) => total_bonds += n,
        None => (),
    }
    let h: i32 = match chemical_formula.get("H") {
        Some(n) => *n,
        None => 0,
    };
    total_bonds = (total_bonds + h) / 2;
    let assigned_bonds = total_bonds - h;
    (total_bonds, assigned_bonds)
}

// PROBABLY DONT NEED THIS
// Computes the index of gen deficiency  to find level of unsaturation
// 1 degree of unsaturation = 1 ring or double bond in final structure
// IHD = (2C + 2 + N - H - X) / 2 where X is halogens
pub fn compute_ihd(elements: &HashMap<&str, i32>) -> i32 {
    let mut ihd: i32 = 2;
    for (k, v) in elements.iter() {
        match k.as_ref() {
            "C" => ihd += 2 * v,
            "H" => ihd -= v,
            "N" => ihd += v,
            "O" => continue,
            "Cl" => ihd -= v,
            "Br" => ihd -= v,
            _ => panic!("Unrecognized element in chemical formula"),
        }
    }
    ihd/2
}
#[test]
fn test_ihd() {
    let caffeine: HashMap<&str, i32> = [("C", 8), ("H", 10),
        ("N", 4), ("O", 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&caffeine), 6);
    let acetic_acid: HashMap<&str, i32> = [("C", 2), ("H", 4),
        ("O", 2)].iter().cloned().collect();
    assert_eq!(compute_ihd(&acetic_acid), 1);
}
