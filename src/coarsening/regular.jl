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

function coarsen(grid::QuasiCartesianGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  RegularGrid(minimum(grid), maximum(grid), dims=size(grid) .÷ factors)
end

function coarsen(grid::RectilinearGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  dims = vsize(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], paramdim(grid))
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> xyzₛ[i][rngs[i]], paramdim(grid))
  RectilinearGrid{manifold(grid),crs(grid)}(xyzₜ)
end

function coarsen(grid::StructuredGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  dims = vsize(grid)
  rngs = ntuple(i -> 1:factors[i]:dims[i], paramdim(grid))
  XYZₛ = XYZ(grid)
  XYZₜ = ntuple(i -> XYZₛ[i][rngs...], paramdim(grid))
  StructuredGrid{manifold(grid),crs(grid)}(XYZₜ)
end

coarsen(grid::TransformedGrid, method::RegularCoarsening) =
  TransformedGrid(coarsen(parent(grid), method), transform(grid))
