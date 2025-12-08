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

function coarsen(grid::RegularGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  dims = _coarsesize(grid, factors)
  orig = minimum(grid)
  cmin = CoordRefSystems.values(coords(minimum(grid)))
  cmax = CoordRefSystems.values(coords(maximum(grid)))
  spac = (cmax .- cmin) ./ dims
  topo = GridTopology(dims, isperiodic(grid))
  RegularGrid(orig, spac, topo)
end

function coarsen(grid::RectilinearGrid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  inds = _coarseinds(grid, factors)
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> xyzₛ[i][inds[i]], paramdim(grid))
  dims = length.(xyzₜ) .- .!isperiodic(grid)
  topo = GridTopology(dims, isperiodic(grid))
  RectilinearGrid{manifold(grid),crs(grid)}(xyzₜ, topo)
end

function coarsen(grid::Grid, method::RegularCoarsening)
  factors = fitdims(method.factors, paramdim(grid))
  inds = _coarseinds(grid, factors)
  XYZₛ = XYZ(grid)
  XYZₜ = ntuple(i -> XYZₛ[i][inds...], paramdim(grid))
  dims = size(first(XYZₜ)) .- .!isperiodic(grid)
  topo = GridTopology(dims, isperiodic(grid))
  StructuredGrid{manifold(grid),crs(grid)}(XYZₜ, topo)
end

coarsen(grid::TransformedGrid, method::RegularCoarsening) =
  TransformedGrid(coarsen(parent(grid), method), transform(grid))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _coarsesize(grid, factors)
  dims = size(grid)
  axes = ntuple(i -> 1:dims[i], paramdim(grid))
  size(TileIterator(axes, factors))
end

_coarsevsize(grid, factors) = _coarsesize(grid, factors) .+ .!isperiodic(grid)

function _coarseinds(grid, factors)
  vdims = vsize(grid)
  cdims = _coarsevsize(grid, factors)
  ntuple(paramdim(grid)) do i
    floor.(Int, range(start=1, stop=vdims[i], length=cdims[i]))
  end
end
