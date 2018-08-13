# Molecule
A C13 NMR interpreter written in Rust. From an input file containing a chemical formula and chemical shifts a genetic algorithm is used to find the interpreted structure.

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

File Input is expected to be two lines. The first containing the chemical formula and the second csv floats for each peak in the C13 spectra including symmetrical peaks.
```
C2H4O2 
162.0, 51.0
```

Currently chemical formulas must explicitly list the number for each element, C2H4O
will not recognize the lone Oxygen.
Elements are limited to C, H, O, N, Cl, Br

Unit tests are run with `Cargo test`

## Chemical Shift Calculation
For a molecule chain ...C&#949;-C&#948;-C&#947;-C&#946;-C&#945;-**C**-C&#945;-C&#946;-C&#947;C&#948;-C&#949;... \
The chemical shift of **C** &#948;C = -2.3 + 9.1&#945; + 9.4&#946; - 2.5&#947; + 0.3&#948; + 0.1&#949; + &#931; (steric corrections) ppm

### Steric corrections
| Carbon Atom Observed | Primary | Secondary | Tertiary | Quaternary |
|----------------------|:-------:|:---------:|:--------:|-----------:|
| Primary              | 0       | 0         | -1.1     | -3.4       |
| Secondary            | 0       | 0         | -2.5     | -7.5       |
| Tertiary             | 0       | -3.7      | -8.5     | -10.0      |
| Quaternary           | 0       | -8.4      | -10.0    | -12.5      |

&#948;C is further affected by substituent effects too numerous to list here. \
Reference: Introduction to Spectroscopy 4th ed, Pavia, Lampman, Kriz, Vyvyan. Appendix 8.



This program is licensed under the "MIT License". Please see the file LICENSE in
the source distribution of this software for license terms.

Copyright (c) 2018 Jeff Lund
