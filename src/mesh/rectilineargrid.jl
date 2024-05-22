# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RectilinearGrid(x, y, z, ...)

A rectilinear grid with vertices at sorted coordinates `x`, `y`, `z`, ...

## Examples

Create a 2D rectilinear grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> x = 0.0:0.2:1.0
julia> y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
julia> RectilinearGrid(x, y)
```
"""
struct RectilinearGrid{Dim,V<:AbstractVector{<:Len}} <: Grid{Dim}
  xyz::NTuple{Dim,V}
  topology::GridTopology{Dim}
  RectilinearGrid{Dim,V}(xyz, topology) where {Dim,V<:AbstractVector{<:Len}} = new(xyz, topology)
end

function RectilinearGrid(xyz::NTuple{Dim,V}, topology::GridTopology{Dim}) where {Dim,V<:AbstractVector{<:Len}}
  coords = float.(xyz)
  RectilinearGrid{Dim,eltype(coords)}(coords, topology)
end

RectilinearGrid(xyz::NTuple{Dim,V}, topology::GridTopology{Dim}) where {Dim,V<:AbstractVector} =
  RectilinearGrid(addunit.(xyz, u"m"), topology)

function RectilinearGrid(xyz::Tuple)
  coords = promote(collect.(xyz)...)
  topology = GridTopology(length.(coords) .- 1)
  RectilinearGrid(coords, topology)
end

RectilinearGrid(xyz...) = RectilinearGrid(xyz)

lentype(::Type{<:RectilinearGrid{Dim,V}}) where {Dim,V} = eltype(V)

vertex(g::RectilinearGrid{Dim}, ijk::Dims{Dim}) where {Dim} = Point(getindex.(g.xyz, ijk))

xyz(g::RectilinearGrid) = g.xyz

XYZ(g::RectilinearGrid) = XYZ(xyz(g))

function centroid(g::RectilinearGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = vertex(g, ijk)
  p2 = vertex(g, ijk .+ 1)
  Point(coords((to(p1) + to(p2)) / 2))
end

function Base.getindex(g::RectilinearGrid{Dim}, I::CartesianIndices{Dim}) where {Dim}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  start = Tuple(first(I))
  stop = Tuple(last(I)) .+ 1
  xyz = ntuple(i -> g.xyz[i][start[i]:stop[i]], Dim)
  RectilinearGrid(xyz, GridTopology(dims))
end

function Base.summary(io::IO, g::RectilinearGrid)
  join(io, size(g), "×")
  print(io, " RectilinearGrid")
end
