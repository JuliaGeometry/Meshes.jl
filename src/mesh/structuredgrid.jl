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
struct StructuredGrid{Datum,Dim,A<:AbstractArray{<:Len}} <: Grid{Dim}
  XYZ::NTuple{Dim,A}
  topology::GridTopology{Dim}
  StructuredGrid{Datum,Dim,A}(XYZ, topology) where {Datum,Dim,A<:AbstractArray{<:Len}} = new(XYZ, topology)
end

function StructuredGrid{Datum}(
  XYZ::NTuple{Dim,A},
  topology::GridTopology{Dim}
) where {Datum,Dim,A<:AbstractArray{<:Len}}
  coords = float.(XYZ)
  StructuredGrid{Datum,Dim,eltype(coords)}(coords, topology)
end

StructuredGrid{Datum}(XYZ::NTuple{Dim,A}, topology::GridTopology{Dim}) where {Datum,Dim,A<:AbstractArray} =
  StructuredGrid{Datum}(addunit.(XYZ, u"m"), topology)

function StructuredGrid{Datum}(XYZ::Tuple) where {Datum}
  coords = promote(XYZ...)
  topology = GridTopology(size(first(coords)) .- 1)
  StructuredGrid{Datum}(coords, topology)
end

StructuredGrid{Datum}(XYZ...) where {Datum} = StructuredGrid{Datum}(XYZ)

StructuredGrid(args...) = StructuredGrid{NoDatum}(args...)

lentype(::Type{<:StructuredGrid{Datum,Dim,A}}) where {Datum,Dim,A} = eltype(A)

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
