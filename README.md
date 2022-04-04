[//]: Logo
<p align="center">
<img
    src="./docs/src/assets/logo256px.svg"
    width=256px
    >
</p>

# Integrated Residuals Solver
[//]: Badges
[![Stable](https://img.shields.io/badge/docs-v0.1-blue.svg)](https://judo-dev.github.io/Interesso.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://judo-dev.github.io/Interesso.jl/dev)
[![Build Status](https://github.com/JuDO-dev/Interesso.jl/workflows/CI/badge.svg)](https://github.com/JuDO-dev/Interesso.jl/actions)
[![Coverage](https://codecov.io/gh/JuDO-dev/Interesso.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuDO-dev/Interesso.jl)

[//]: Description
Solves **Dynamic Feasibility** and **Dynamic Optimisation** problems (think Optimal Control) using an Integrated Residuals method.
## Installation

```julia
julia> ]
  pkg> add Interesso
```

## Example using solve()
```julia
using Interesso

dop = DOProblem();

solution = solve(dop, LeastSquares());
```

## Example using iterator()
```julia
using Interesso

dop = DOProblem();

I = DOIterator(dop, LeastSquares());

collect(I);
```