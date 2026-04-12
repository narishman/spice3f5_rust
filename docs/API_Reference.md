# API Reference Information

`spice-rs` is principally a mathematical numerical framework natively designed to allow flexible usage. Users can bypass the text parser fully using the structural APIs intrinsically generating equivalent parameter boundaries, or simply utilize the parser independently for programmatic textual evaluations natively.

## Circuit Object Instantiation

The entire logic sequence revolves heavily around an intrinsic localized standard variable: `Circuit`. The circuit evaluates the actual N-R sequences directly internally matching boundaries structurally natively.

```rust
use spice_core::circuit::Circuit;
let mut circuit = Circuit::new();

// Dynamically assign or locate mathematical variable indexes safely maps to the grid boundary recursively natively
let node_a = circuit.get_or_create_node("NodeIdentifier_A", spice_core::node::NodeType::Voltage);
let node_gnd = 0; // Evaluates inherently mathematically to structural zeroes across mappings
```

## Adding Component Boundaries

Structural components actively exist inside a generic `enum` wrapping individual memory states: `DeviceEnum`.

```rust
use spice_core::devices::{DeviceEnum, resistor::Resistor, mos1::{Mos1, MosType}};
let mut devices = vec![
    DeviceEnum::Resistor(Resistor::new("R1", node_a, node_gnd, 1000.0)),
    // Dynamically insert Complex boundaries tracking iterative equivalents natively!
    DeviceEnum::Mos1(Mos1::new("M1_Fet", MosType::NMOS, drain, gate, source, bulk))
];
```

## Parsing Engine Integrations

The standalone AST parsing boundaries uniquely trace string operations into mathematical definitions cleanly avoiding unsafe generic boundaries identically passing strict Rust syntax configurations naturally:

```rust
use spice_parse::{lexer::Lexer, parser::Parser};
use spice_core::loader::load_from_ast;

let internal_ast_token_string = "
M1 1 2 0 0
R1 1 0 1000.0
";

let tokens = Lexer::new(internal_ast_token_string).tokenize();
let statements = Parser::new(tokens).parse().unwrap();

// Wire bindings perfectly across dynamic parameters automatically translating components recursively structurally!
let active_components = load_from_ast(&mut circuit, statements).unwrap();
```

## Operational Evaluators

Evaluations natively run recursive loops mapping output mathematical variable matrices natively calculating outputs!

```rust
use spice_core::ni::do_dc_op;

// Execute analytical limits boundaries inherently resolving structures natively accurately
do_dc_op(&mut circuit, &mut devices).expect("N-R Integration mathematical bounds failed natively!");
```
