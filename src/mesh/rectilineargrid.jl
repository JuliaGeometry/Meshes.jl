# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RectilinearGrid(x, y, z, ...)
    RectilinearGrid{Datum}(x, y, z, ...)

A rectilinear grid with vertices at sorted coordinates `x`, `y`, `z`, ...,
and a given `Datum` (default to `NoDatum`).

## Examples

Create a 2D rectilinear grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> x = 0.0:0.2:1.0
julia> y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
julia> RectilinearGrid(x, y)
```
"""
struct RectilinearGrid{Datum,Dim,ℒ<:Len,V<:AbstractVector{ℒ}} <: Grid{Dim,Cartesian{Datum,Dim,ℒ}}
  xyz::NTuple{Dim,V}
  topology::GridTopology{Dim}

  function RectilinearGrid{Datum}(
    xyz::NTuple{Dim,<:AbstractVector{<:Len}},
    topology::GridTopology{Dim}
  ) where {Datum,Dim}
    coords = float.(xyz)
    V = eltype(coords)
    new{Datum,Dim,eltype(V),V}(coords, topology)
  end
end

RectilinearGrid{Datum}(xyz::NTuple{Dim,<:AbstractVector}, topology::GridTopology{Dim}) where {Datum,Dim} =
  RectilinearGrid{Datum}(addunit.(xyz, u"m"), topology)

function RectilinearGrid{Datum}(xyz::Tuple) where {Datum}
  coords = promote(collect.(xyz)...)
  topology = GridTopology(length.(coords) .- 1)
  RectilinearGrid{Datum}(coords, topology)
end

RectilinearGrid{Datum}(xyz...) where {Datum} = RectilinearGrid{Datum}(xyz)

RectilinearGrid(args...) = RectilinearGrid{NoDatum}(args...)

vertex(g::RectilinearGrid{Datum,Dim}, ijk::Dims{Dim}) where {Datum,Dim} = Point(Cartesian{Datum}(getindex.(g.xyz, ijk)))

xyz(g::RectilinearGrid) = g.xyz

XYZ(g::RectilinearGrid) = XYZ(xyz(g))

function centroid(g::RectilinearGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = vertex(g, ijk)
  p2 = vertex(g, ijk .+ 1)
  withdatum(g, (to(p1) + to(p2)) / 2)
end

function Base.getindex(g::RectilinearGrid{Datum,Dim}, I::CartesianIndices{Dim}) where {Datum,Dim}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  start = Tuple(first(I))
  stop = Tuple(last(I)) .+ 1
  xyz = ntuple(i -> g.xyz[i][start[i]:stop[i]], Dim)
  RectilinearGrid{Datum}(xyz, GridTopology(dims))
end

function Base.summary(io::IO, g::RectilinearGrid)
  join(io, size(g), "×")
  print(io, " RectilinearGrid")
end
