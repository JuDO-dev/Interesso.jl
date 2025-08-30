```@meta
CurrentModule = Interesso
```

# Interesso

A powerful and extendable solver for nonlinear **Trajectory Optimization** (aka Optimal Control) problems using pseudo-spectral methods.

## Features
- **Flexible Intervals**: Supports flexible discretization intervals, avoiding mesh *h*-refinement.
- **Tight Constraints**: Uses Bernstein polynomial coefficients to constraint dynamic variables tightly.
- **JuMP Backend**: Choose your favorite optimizer from JuMP.jl's [supported solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

## Installation

```julia
julia> ]
  pkg> add Interesso
```

## License
Interesso.jl is licensed under the MIT License.

## Contributers
 - Eduardo Vila [@e-duar-do](https://github.com/e-duar-do)

We welcome contributions! Create an issue or pull request on our [GitHub repo](https://github.com/JuDO-dev/Interesso.jl).