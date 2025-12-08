# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularRefinement(f₁, f₂, ..., fₙ)

Refine each dimension of the grid by given factors `f₁`, `f₂`, ..., `fₙ`.

## Examples

```julia
refine(grid2D, RegularRefinement(2, 3))
refine(grid3D, RegularRefinement(2, 3, 1))
```
"""
struct RegularRefinement{N} <: RefinementMethod
  factors::Dims{N}
end

RegularRefinement(factors::Vararg{Int,N}) where {N} = RegularRefinement(factors)

function refine(grid::RegularGrid, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  dims = size(grid) .* factors
  orig = minimum(grid)
  cmin = CoordRefSystems.values(coords(minimum(grid)))
  cmax = CoordRefSystems.values(coords(maximum(grid)))
  spac = (cmax .- cmin) ./ dims
  topo = GridTopology(dims, isperiodic(grid))
  RegularGrid(orig, spac, topo)
end

function refine(grid::RectilinearGrid, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  xyzₛ = xyz(grid)
  xyzₜ = ntuple(i -> _refinedims(xyzₛ[i], factors[i]), paramdim(grid))
  dims = length.(xyzₜ) .- 1
  topo = GridTopology(dims, isperiodic(grid))
  RectilinearGrid{manifold(grid),crs(grid)}(xyzₜ, topo)
end

function refine(grid::Grid, method::RegularRefinement)
  factors = fitdims(method.factors, paramdim(grid))
  XYZₜ = XYZ(grid, factors)
  dims = size(first(XYZₜ)) .- 1
  topo = GridTopology(dims, isperiodic(grid))
  StructuredGrid{manifold(grid),crs(grid)}(XYZₜ, topo)
end

refine(grid::TransformedGrid, method::RegularRefinement) =
  TransformedGrid(refine(parent(grid), method), transform(grid))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _refinedims(x, f)
  x′ = mapreduce(vcat, 1:(length(x) - 1)) do i
    range(x[i], x[i + 1], f + 1)[begin:(end - 1)]
  end
  push!(x′, last(x))
  x′
end
