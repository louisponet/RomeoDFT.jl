# States
```@meta
CurrentModule = RomeoDFT
```
A [`State`](@ref) represents a set of occupied local manifolds, one for each ion to which occupation constraints are applied, or that have a Hubbard U parameter.
Usually no direct interaction with these is needed as the search is happening, but they are used for some postprocessing analysis.

They hold:
- `occupations`: the occupation matrices for each ion
- `totoccs`: total electron occupation of each ion
- `magmoms`: magnetic moments defined as `tr(up) - tr(down)`
- `eigvals`: the eigen values for each occupation matrix
- `eigvecs`: the eigen vectors for each occupation matrix

## Euclidean distance

The `Euclidean` distance metric is defined between states, which essentially does
```math
\sqrt{\sum_I \sum_{\alpha\beta}(n^I_{1, \alpha\beta} - n^I_{2, \alpha\beta})^2}
```

and can be called like
```julia
RomeoDFT.Euclidean()(state1, state2)
```

## Refeference
```@docs
State
generate_Hubbard_occupations
generate_starting_ns_eigenvalue
```

