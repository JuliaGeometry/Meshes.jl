# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StructuredGrid(X, Y, Z, ...)
    StructuredGrid{Datum}(X, Y, Z, ...)

A structured grid with vertices at sorted coordinates `X`, `Y`, `Z`, ...,
and a given `Datum` (default to `NoDatum`).

## Examples

Create a 2D structured grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> X = repeat(0.0:0.2:1.0, 1, 6)
julia> Y = repeat([0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
julia> StructuredGrid(X, Y)
```
"""
struct StructuredGrid{Datum,Dim,ℒ<:Len,A<:AbstractArray{ℒ}} <: Grid{Dim,Cartesian{Datum,Dim,ℒ}}
  XYZ::NTuple{Dim,A}
  topology::GridTopology{Dim}
  StructuredGrid{Datum,Dim,ℒ,A}(XYZ, topology) where {Datum,Dim,ℒ<:Len,A<:AbstractArray{ℒ}} = new(XYZ, topology)
end

function StructuredGrid{Datum}(XYZ::NTuple{Dim,<:AbstractArray{<:Len}}, topology::GridTopology{Dim}) where {Datum,Dim}
  coords = float.(XYZ)
  A = eltype(coords)
  StructuredGrid{Datum,Dim,eltype(A),A}(coords, topology)
end

StructuredGrid{Datum}(XYZ::NTuple{Dim,<:AbstractArray}, topology::GridTopology{Dim}) where {Datum,Dim} =
  StructuredGrid{Datum}(addunit.(XYZ, u"m"), topology)

function StructuredGrid{Datum}(XYZ::Tuple) where {Datum}
  coords = promote(XYZ...)
  topology = GridTopology(size(first(coords)) .- 1)
  StructuredGrid{Datum}(coords, topology)
end

StructuredGrid{Datum}(XYZ...) where {Datum} = StructuredGrid{Datum}(XYZ)

StructuredGrid(args...) = StructuredGrid{NoDatum}(args...)

vertex(g::StructuredGrid{Datum,Dim}, ijk::Dims{Dim}) where {Datum,Dim} =
  Point(Cartesian{Datum}(ntuple(d -> g.XYZ[d][ijk...], Dim)))

XYZ(g::StructuredGrid) = g.XYZ

function Base.getindex(g::StructuredGrid{Datum,Dim}, I::CartesianIndices{Dim}) where {Datum,Dim}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
  XYZ = ntuple(i -> g.XYZ[i][cinds], Dim)
  StructuredGrid{Datum}(XYZ, GridTopology(dims))
end

function Base.summary(io::IO, g::StructuredGrid)
  join(io, size(g), "×")
  print(io, " StructuredGrid")
end
