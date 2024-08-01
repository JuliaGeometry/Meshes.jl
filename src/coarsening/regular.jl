# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularCoarsening(f₁, f₂, ..., fₙ)

Coarsen each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

## Examples

```julia
coarsen(grid2D, RegularCoarsening(2, 3))
coarsen(grid3D, RegularCoarsening(2, 3, 1))
```
"""
struct RegularCoarsening{N} <: CoarseningMethod
  factors::Dims{N}
end

RegularCoarsening(factors::Vararg{Int,N}) where {N} = RegularCoarsening(factors)

function coarsen(grid::CartesianGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, embeddim(grid))
  CartesianGrid(minimum(grid), maximum(grid), dims=size(grid) .÷ factors)
end

function coarsen(grid::RectilinearGrid{M,C}, method::RegularCoarsening) where {M,C}
  factors = fitdims(method.factors, embeddim(grid))
  dims = size(grid) .+ .!isperiodic(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], embeddim(grid))
  xyzₛ = map(x -> ustrip.(x), xyz(grid))
  xyzₜ = ntuple(i -> xyzₛ[i][rngs[i]], embeddim(grid))
  RectilinearGrid{M,C}(xyzₜ)
end

function coarsen(grid::StructuredGrid{Datum}, method::RegularCoarsening) where {Datum}
  factors = fitdims(method.factors, embeddim(grid))
  dims = size(grid) .+ .!isperiodic(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], embeddim(grid))
  XYZₛ = XYZ(grid)
  XYZₜ = ntuple(i -> XYZₛ[i][rngs...], embeddim(grid))
  StructuredGrid{Datum}(XYZₜ)
end
