# spice-rs Wiki

Welcome to the central documentation hub for `spice-rs`. Here, we provide an extensive deep dive into the underlying architecture, data structures, and mathematical porting strategies utilized to perfectly recreate the historical Berkeley SPICE 3 in a fully modern, memory-safe Rust codebase.

## Documentation Index

- [Architecture Overview](Architecture.md) - A structural deep-dive into the workspace crates, simulation lifecycle boundaries, and code responsibilities.
- [Mathematical Engine (MNA & Sparse LU)](Mathematical_Engine.md) - A rigorous walkthrough of Modified Nodal Analysis (MNA), how generic structures map into the computational space, and exactly how the specialized Sparse Matrix handles factorization and memory.
- [Non-Linear Active Mechanics](NonLinear_Active_Mechanics.md) - Explanation of the Newton-Raphson (N-R) dynamic solver bounding rules, non-linear iterative evaluation limits, and simulation convergence mathematics.
- [MOSFET Dev-Model Deep-Dive](MOSFET_Physics.md) - Specialized exploration evaluating the Shichman-Hodges Level-1 MOSFET (and non-linear component derivation logic mapping) down to the algorithm. 
- [API Reference Guide](API_Reference.md) - Practical programmatic developer information on manually injecting components and binding the netlist ASTs.
