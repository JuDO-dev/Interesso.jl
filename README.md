[//]: Logo
<p align="center">
<img
    src="./docs/src/assets/logo.svg"
    width=128px
    >
</p>

# Interesso
[//]: Badges
[![Stable](https://img.shields.io/badge/docs-v0.1-blue.svg)](https://judo-dev.github.io/Interesso.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://judo-dev.github.io/Interesso.jl/dev)
[![Build Status](https://github.com/JuDO-dev/Interesso.jl/workflows/CI/badge.svg)](https://github.com/JuDO-dev/Interesso.jl/actions)
[![Coverage](https://codecov.io/gh/JuDO-dev/Interesso.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuDO-dev/Interesso.jl)

[//]: Description
A powerful and extendable solver for nonlinear **Trajectory Optimization** (aka Optimal Control) problems using pseudo-spectral methods.

## Highlights
- **Flexible Intervals**: Supports flexible discretization intervals, avoiding mesh *h*-refinement.
- **Tight Constraints**: Uses Bernstein polynomial coefficients to constraint dynamic variables tightly.
- **JuMP Backend**: Choose your favorite optimizer from JuMP.jl's [supported solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

## Installation

```julia
julia> ]
  pkg> add Interesso
```

## Usage
Find usage examples in the [documentation](https://judo-dev.github.io/Interesso.jl).