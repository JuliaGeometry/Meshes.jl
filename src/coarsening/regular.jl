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

function coarsen(grid::OrthoRegularGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  dims = _targetsize(grid, factors)
  RegularGrid(minimum(grid), maximum(grid), dims=dims)
end

function coarsen(grid::RectilinearGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  inds = _targetinds(grid, factors)
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> xyzₛ[i][inds[i]], paramdim(grid))
  RectilinearGrid{manifold(grid),crs(grid)}(xyzₜ)
end

function coarsen(grid::StructuredGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  inds = _targetinds(grid, factors)
  XYZₛ = XYZ(grid)
  XYZₜ = ntuple(i -> XYZₛ[i][inds...], paramdim(grid))
  StructuredGrid{manifold(grid),crs(grid)}(XYZₜ)
end

coarsen(grid::TransformedGrid, method::RegularCoarsening) =
  TransformedGrid(coarsen(parent(grid), method), transform(grid))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _targetsize(grid, factors)
  dims = size(grid)
  axes = ntuple(i -> 1:dims[i], paramdim(grid))
  size(TileIterator(axes, factors))
end

_targetvsize(grid, factors) = _targetsize(grid, factors) .+ .!isperiodic(grid)

function _targetinds(grid, factors)
  dims = vsize(grid)
  tdims = _targetvsize(grid, factors)
  ntuple(i -> floor.(Int, range(start=1, stop=dims[i], length=tdims[i])), paramdim(grid))
end
