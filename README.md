# spice-rs

> **Note:** This project is a direct 1:1 algorithmic port of the original Berkeley SPICE 3f5, accomplished using [Google Antigravity](https://deepmind.google/) with very little human oversight. 

`spice-rs` is a high-performance, 100% memory-safe Rust modernization of the classic Berkeley SPICE 3f5 circuit simulator. 

This project aims to deliver the exact same numerical precision and analytical robustness of SPICE 3 in a completely safe, modular, and performant Rust architecture. It is built as a complete 1:1 algorithmic port of the mathematical engine, sparse matrix solver, and non-linear component evaluation pipelines found in SPICE.

## Features

- **Memory-Safe Sparse LU Solver**: Identical 1:1 port of Kenneth Kundert’s `SPARSE` 1.3 linear equation solver, completely modernized to avoid unsafe pointer arithmetic using safe, contiguous Arena allocations (`spice_math`).
- **Dynamic Newton-Raphson Iteration**: Fully implemented DC Operating Point (`do_dc_op`) and Transient (`do_tran`) integration loops matching original SPICE limits and tolerance bounding (`spice_core::ni`).
- **Active Semiconductor Logic**: 
    - Ebers-Moll / Gummel-Poon BJT (`Bjt`)
    - Shichman-Hodges Level-1 MOSFET (`Mos1`)
    - PN Junction Diode (`Diode`)
- **Passive & Linear Netlist Models**: Fully integrated `Resistor`, `Capacitor`, `Inductor`, `VoltageSource`, and `CurrentSource`.
- **AST Parsing Engine**: Reads standard SPICE netlists into executable AST nodes via `spice_parse`.

## Project Structure

The workspace is divided into three distinct crates to ensure clean separation of concerns:

- `spice_math`: A standalone, highly tailored sparse matrix solver heavily optimized for nodal admittance matrices. Calculates Markowitz factorization and solves the matrix.
- `spice_parse`: The lexer and parser engine responsible for breaking text netlists down into ASTs.
- `spice_core`: The main interface logic. It instantiates the `Circuit`, manages integration history states, dynamically assigns evaluation limits across non-linear junctions, and loads the matrix via specific device structural equations.

## Documentation Wiki

For a deep dive into the underlying architecture, data structures, and mathematical porting strategies utilized to perfectly recreate the historical Berkeley SPICE 3, please check out the **[spice-rs Wiki](docs/Home.md)**!

The Wiki contains complete detailed explanations of:
- Architectural flow of components.
- Modified Nodal Analysis and Sparse Matrix Formulations heavily diagrammed.
- Newton-Raphson mathematics bounding active iterative structures.
- A highly detailed Shichman-Hodges Level 1 MOSFET algorithmic sequence exploration.
- Robust programmatic API tutorials.

## Installation

Ensure you have Rust installed via `rustup`. Since `spice-rs` is a standard Cargo workspace, installation is trivial:

```bash
git clone https://github.com/yourusername/spice-rs.git
cd spice-rs
cargo build --release
```

## Basic Usage

The engine is primarily designed to be imported as a library for high-speed simulation backends. Below is an example of creating a simple nodal circuit network mathematically inline.

### Constructing a Basic Circuit

```rust
use spice_core::circuit::Circuit;
use spice_core::devices::{DeviceEnum, resistor::Resistor, vsource::VoltageSource};
use spice_core::ni::do_dc_op;

fn main() {
    let mut circuit = Circuit::new();
    
    // In SPICE, Node 0 is always Ground
    let gnd = 0;
    
    let n1 = circuit.get_or_create_node("1", spice_core::node::NodeType::Voltage);
    let n2 = circuit.get_or_create_node("2", spice_core::node::NodeType::Voltage);

    let mut devices = vec![
        DeviceEnum::VoltageSource(VoltageSource::new("V1", n1, gnd, 5.0)),
        DeviceEnum::Resistor(Resistor::new("R1", n1, n2, 1000.0)),
        DeviceEnum::Resistor(Resistor::new("R2", n2, gnd, 1000.0)),
    ];
    
    for dev in devices.iter_mut() {
        dev.setup(&mut circuit).unwrap();
    }

    // Solve for the DC Operating Point analytically
    do_dc_op(&mut circuit, &mut devices).expect("NR Convergence Failed");

    println!("Node 1 Voltage: {} V", circuit.rhs[n1]); // 5.0V
    println!("Node 2 Voltage: {} V", circuit.rhs[n2]); // 2.5V (Voltage Divider)
}
```

### Loading from an AST
Alternatively, if you want to bypass manual structuring and load from a parsed string:

```rust
use spice_core::circuit::Circuit;
use spice_parse::{lexer::Lexer, parser::Parser};
use spice_core::loader::load_from_ast;

fn main() {
    let input = "
    V1 1 0 5.0
    R1 1 2 1000.0
    R2 2 0 1000.0
    ";
    
    let mut lexer = Lexer::new(input);
    let tokens = lexer.tokenize();
    let mut parser = Parser::new(tokens);
    let ast_stmts = parser.parse().unwrap();
    
    let mut circuit = Circuit::new();
    let devices = load_from_ast(&mut circuit, ast_stmts).unwrap();
    
    // Proceed to solve with `do_dc_op` or `do_tran`
}
```

## Running Tests

All numerical logic has been integrated with pure physics proofs validating adherence to original SPICE matrix derivations. You can run these internal validations at any point:

```bash
cargo test --workspace
```

## Contributing
As `spice-rs` is a strict math port, contributing components and models should mimic the parameterization and limit bounds implemented by Berkeley. All contributions focusing on porting level 2+ MOSFETs (e.g., BSIM3/4) must conform to standard Sparse-LU topological loading strategies.

## License
`spice-rs` operates under the identical BSD Copyright as the original Berkeley SPICE 3 source code it was derived from. 

Copyright (c) 1985-1991 The Regents of the University of California.
All rights reserved.

See the `LICENSE` file for more details.
