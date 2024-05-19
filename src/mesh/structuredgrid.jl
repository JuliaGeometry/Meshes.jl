# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StructuredGrid(X, Y, Z, ...)

A structured grid with vertices at sorted coordinates `X`, `Y`, `Z`, ...

## Examples

Create a 2D structured grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> X = repeat(0.0:0.2:1.0, 1, 6)
julia> Y = repeat([0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
julia> StructuredGrid(X, Y)
```
"""
struct StructuredGrid{Dim,A<:AbstractArray{<:Len}} <: Grid{Dim}
  XYZ::NTuple{Dim,A}
  topology::GridTopology{Dim}
  StructuredGrid{Dim,A}(XYZ, topology) where {Dim,A<:AbstractArray{<:Len}} = new(XYZ, topology)
end

function StructuredGrid(XYZ::NTuple{Dim,A}, topology::GridTopology{Dim}) where {Dim,A<:AbstractArray{<:Len}}
  coords = float.(XYZ)
  StructuredGrid{Dim,eltype(coords)}(coords, topology)
end

StructuredGrid(XYZ::NTuple{Dim,A}, topology::GridTopology{Dim}) where {Dim,A<:AbstractArray} =
  StructuredGrid(addunit.(XYZ, u"m"), topology)

function StructuredGrid(XYZ::Tuple)
  coords = promote(XYZ...)
  topology = GridTopology(size(first(coords)) .- 1)
  StructuredGrid(coords, topology)
end

StructuredGrid(XYZ...) = StructuredGrid(XYZ)

lentype(::Type{<:StructuredGrid{Dim,A}}) where {Dim,A} = eltype(A)

vertex(g::StructuredGrid{Dim}, ijk::Dims{Dim}) where {Dim} = Point(ntuple(d -> g.XYZ[d][ijk...], Dim))

XYZ(g::StructuredGrid) = g.XYZ

function Base.getindex(g::StructuredGrid{Dim}, I::CartesianIndices{Dim}) where {Dim}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
  XYZ = ntuple(i -> g.XYZ[i][cinds], Dim)
  StructuredGrid(XYZ, GridTopology(dims))
end

function Base.summary(io::IO, g::StructuredGrid)
  join(io, size(g), "×")
  print(io, " StructuredGrid")
end
