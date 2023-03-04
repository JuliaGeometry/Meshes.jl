# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RectilinearGrid(x, y, z, ...)

A rectilinear grid with vertices at coordinates `x`, `y`, `z`, ...

## Examples

Create a 2D rectilinear grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> x = 0.0:0.2:1.0
julia> y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
julia> RectilinearGrid(x, y)
```
"""
struct RectilinearGrid{Dim,T,V<:AbstractVector{T}} <: Grid{Dim,T}
  xyz::NTuple{Dim,V}
  topology::GridTopology{Dim}
end

function RectilinearGrid(xyz::Tuple)
  coords = promote(collect.(xyz)...)
  topology = GridTopology(length.(coords) .- 1)
  return RectilinearGrid(coords, topology)
end

RectilinearGrid(xyz...) = RectilinearGrid(xyz)

cart2vert(g::RectilinearGrid, ijk::Tuple) = Point(getindex.(g.xyz, ijk))

function centroid(g::RectilinearGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = cart2vert(g, ijk)
  p2 = cart2vert(g, ijk .+ 1)
  return Point((coordinates(p1) + coordinates(p2)) / 2)
end
