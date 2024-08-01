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
struct RectilinearGrid{M<:Manifold,C<:CRS,Dim,V<:AbstractVector} <: Grid{M,C,Dim}
  xyz::NTuple{Dim,V}
  topology::GridTopology{Dim}

  function RectilinearGrid{M,C}(xyz::NTuple{Dim,<:AbstractVector}, topology::GridTopology{Dim}) where {M,C,Dim}
    coords = float.(xyz)
    V = eltype(coords)
    new{M,C,Dim,eltype(V),V}(coords, topology)
  end
end

function RectilinearGrid{M,C}(xyz::Tuple) where {M,C}
  coords = promote(collect.(xyz)...)
  topology = GridTopology(length.(coords) .- 1)
  RectilinearGrid{M,C}(coords, topology)
end

RectilinearGrid{M,C}(xyz...) where {M,C} = RectilinearGrid{M,C}(xyz)

function RectilinearGrid(xyz::Tuple)
  Dim = length(xyz)
  C = Cartesian{NoDatum,Dim,Met{Float64}}
  RectilinearGrid{𝔼{Dim},C}(args...)
end

RectilinearGrid(xyz...) = RectilinearGrid(xyz)

vertex(g::RectilinearGrid{M,C}, ijk::Dims) where {M,C} = Point(CoordRefSystems.reconstruct(C, getindex.(g.xyz, ijk)))

xyz(g::RectilinearGrid{M,C}) where {M,C} = g.xyz .* CoordRefSystems.units(C)

XYZ(g::RectilinearGrid) = XYZ(xyz(g))

function centroid(g::RectilinearGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = vertex(g, ijk)
  p2 = vertex(g, ijk .+ 1)
  withcrs(g, (to(p1) + to(p2)) / 2)
end

function Base.getindex(g::RectilinearGrid{M,C}, I::CartesianIndices) where {M,C}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  start = Tuple(first(I))
  stop = Tuple(last(I)) .+ 1
  xyz = ntuple(i -> g.xyz[i][start[i]:stop[i]], embeddim(g))
  RectilinearGrid{M,C}(xyz, GridTopology(dims))
end

function Base.summary(io::IO, g::RectilinearGrid)
  join(io, size(g), "×")
  print(io, " RectilinearGrid")
end
