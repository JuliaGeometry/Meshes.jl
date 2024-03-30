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

function coarsen(grid::CartesianGrid{Dim}, method::RegularCoarsening) where {Dim}
  factors = fitdims(method.factors, Dim)
  CartesianGrid(minimum(grid), maximum(grid), dims=size(grid) .÷ factors)
end

function coarsen(grid::RectilinearGrid{Dim}, method::RegularCoarsening) where {Dim}
  factors = fitdims(method.factors, Dim)
  dims = size(grid) .+ .!isperiodic(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], Dim)
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> xyzₛ[i][rngs[i]], Dim)
  RectilinearGrid(xyzₜ)
end

function coarsen(grid::StructuredGrid{Dim}, method::RegularCoarsening) where {Dim}
  factors = fitdims(method.factors, Dim)
  dims = size(grid) .+ .!isperiodic(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], Dim)
  XYZₛ = XYZ(grid)
  XYZₜ = ntuple(i -> XYZₛ[i][rngs...], Dim)
  StructuredGrid(XYZₜ)
end
