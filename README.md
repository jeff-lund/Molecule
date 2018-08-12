# Molecule
A C13 NMR interpreter written in Rust. From an input file containing  
a chemical formula and chemical shifts a genetic algorithm is used to
find the interpreted structure.

## Getting Started
Molecule uses the [standard cargo building process](https://doc.rust-lang.org/cargo/guide/working-on-an-existing-project.html).
```sh
git clone https://github.com/jeff-lund/Molecule.git
cd Molecule
cargo build
```

Molecule is run from the command line on a file. Samples files are included in the
`test_files` directory with sub-directories for expected runtimes.
```
cargo run test_files/short/acetic_acid.txt
```

File Input is expected to be two lines.
[Chemical formula i.e C2H4O2] \
[Comma separated chemical shifts i.e 162.0, 51.0] \

Currently chemical formulas must explicitly list the number for each element, C2H4O
will not recognize the lone Oxygen.
Elements are limited to C, H, O, N, Cl, Br

Unit tests are run with `Cargo test`

This program is licensed under the "MIT License". Please see the file LICENSE in
the source distribution of this software for license terms.

Copyright (c) 2018 Jeff Lund
