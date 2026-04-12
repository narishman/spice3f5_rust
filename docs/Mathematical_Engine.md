# The Mathematical Engine

`spice-rs` solves circuits based entirely on Modified Nodal Analysis (MNA), an extension of standard Kirchhoff Current Law matrix combinations that accommodates ideal voltage components and dependent mechanisms inherently.

## MNA Representation

Any valid circuit component can be mapped iteratively mapped down to a generic row and column offset within a multidimensional admittance matrix ($G$):

$$ [G] \cdot [V] = [I] $$

Where $G$ represents admittance vectors across given nodes, and $I$ maps explicit or internally bounded equivalent constant currents feeding the components (such as constant voltage evaluation injections). This results in evaluating matrices mathematically spanning $n \times n$ where $n$ scales identically with the amount of parsed nodes (excluding ideal Ground $0$).

```mermaid
graph LR
    A((Node 1)) -- Resistor --> B((Node 2))
    
    style A fill:#4CAF50,stroke:#388E3C,stroke-width:2px,color:white
    style B fill:#2196F3,stroke:#1976D2,stroke-width:2px,color:white
```

In the diagram above, the connection logically establishes elements in the mathematical grid:
- Positive equivalent conductances sum at $[N_1,N_1]$ and $[N_2,N_2]$.
- Reciprocal negative offsets define $[N_1,N_2]$ and $[N_2,N_1]$.

This exact principle spans dynamically to non-linear components (MOSFETS, Diodes). Every complex branch evaluates standard Jacobian admittances iteratively onto these precise mapping points across the evaluation boundary in `load()` arrays!

## Matrix Structures & Algorithms
### The Memory Problem
Most cells inside a fully realized circuit evaluation matrix are mathematically equivalent to $0$. A 100-node matrix spans 10,000 mapping parameters. If 90% are $0$, statically generating mathematical inversion bounds natively incurs a brutal cubic $O(n^3)$ latency evaluating empty sequences. SPICE circumvented this via linked-list memory structures (`SPARSE`), heavily dynamically traversing parameters at runtime to guarantee minimal compute waste. 

### The Rust Implementation
`spice-rs` completely port-replicates Kenneth Kundert’s initial `SPARSE` 1.3 LU strategies seamlessly without suffering unsafe pointer iterations! 
We utilize a heavily typed, statically structured **Sparse Engine (`spice_math`)**:

- **Arena Based Indexing:** Instead of standard linked-address memory, structures dynamically resolve around statically allocated generic vector indexes avoiding bounds panics safely. Nodes refer to contiguous identifiers natively passing pointer validity boundaries statically across compilation times.

### Markowitz Decomposition
Finding numerical inversions inherently yields $Ax = b$ solutions sequentially via Gaussian LU Elimination formats:
1. Matrix boundaries continuously check mathematical fill-in penalties (`update_markowitz_numbers`).
2. Pivot mapping boundaries calculate optimal element evaluation sequences avoiding mathematical divide-by-zero singularities.
3. Linear boundaries process factorizations cleanly across pre-assigned optimal structural sequences down optimal vectors!
4. Solvers trace linear RHS substitutions.
