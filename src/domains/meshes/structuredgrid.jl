# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StructuredGrid(X, Y, Z, ...)
    StructuredGrid{M,C}(X, Y, Z, ...)

A structured grid with vertices at sorted coordinates `X`, `Y`, `Z`, ...,
manifold `M` (default to `𝔼`) and CRS type `C` (default to `Cartesian`).

    StructuredGrid((X, Y, Z, ...), topology)
    StructuredGrid{M,C}((X, Y, Z, ...), topology)

Alternatively, construct a structured grid with `(X, Y, Z, ...)` coordinates
and grid `topology`. This method is available for advanced use cases involving
periodic dimensions. See [`GridTopology`](@ref) for more details.

## Examples

Create a 2D structured grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> X = repeat(0.0:0.2:1.0, 1, 6)
julia> Y = repeat([0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
julia> StructuredGrid(X, Y)
```
"""
struct StructuredGrid{M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractArray},D,TP<:GridTopology{D}} <: Grid{M,C,D}
  XYZ::X
  topology::TP

  function StructuredGrid{M,C,N,X,D,TP}(
    XYZ,
    topology
  ) where {M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractArray},D,TP<:GridTopology{D}}
    new(XYZ, topology)
  end
end

function StructuredGrid{M,C}(XYZ::NTuple{N,AbstractArray}, topology::GridTopology{D}) where {M<:Manifold,C<:CRS,N,D}
  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)
  nd = ndims(first(XYZ))

  allequal(size(X) for X in XYZ) || throw(ArgumentError("all coordinate arrays must have the same size"))
  N == nc || throw(ArgumentError("CRS with $nc coordinates and $N coordinate arrays provided"))
  D == nd || throw(ArgumentError("$D-dimensional grid topology and $nd-dimensional coordinate arrays provided"))

  XYZ′ = ntuple(i -> numconvert.(T, withunit.(XYZ[i], us[i])), nc)

  StructuredGrid{M,C,N,typeof(XYZ′),D,typeof(topology)}(XYZ′, topology)
end

function StructuredGrid{M,C}(XYZ::NTuple{N,AbstractArray}) where {M<:Manifold,C<:CRS,N}
  topology = GridTopology(size(first(XYZ)) .- 1)
  StructuredGrid{M,C}(XYZ, topology)
end

StructuredGrid{M,C}(XYZ::AbstractArray...) where {M<:Manifold,C<:CRS} = StructuredGrid{M,C}(XYZ)

function StructuredGrid(XYZ::NTuple{N,AbstractArray}, topology::GridTopology) where {N}
  L = promote_type(ntuple(i -> aslentype(float(eltype(XYZ[i]))), N)...)
  M = 𝔼{N}
  C = Cartesian{NoDatum,N,L}
  StructuredGrid{M,C}(XYZ, topology)
end

function StructuredGrid(XYZ::NTuple{N,AbstractArray}) where {N}
  topology = GridTopology(size(first(XYZ)) .- 1)
  StructuredGrid(XYZ, topology)
end

StructuredGrid(XYZ::AbstractArray...) = StructuredGrid(XYZ)

function vertex(g::StructuredGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  Point(ctor(ntuple(d -> g.XYZ[d][ijk...], length(g.XYZ))...))
end

XYZ(g::StructuredGrid) = g.XYZ

@generated function Base.getindex(g::StructuredGrid{M,C,N}, I::CartesianIndices) where {M,C,N}
  exprs = ntuple(i -> :(g.XYZ[$i][cinds]), N)

  quote
    @boundscheck _checkbounds(g, I)
    dims = size(I)
    cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
    XYZ = ($(exprs...),)
    StructuredGrid{M,C}(XYZ, GridTopology(dims))
  end
end

function Base.summary(io::IO, g::StructuredGrid)
  join(io, size(g), "×")
  print(io, " StructuredGrid")
end

function _assertXYZ(XYZ)
  if !allequal(size(X) for X in XYZ)
    throw(ArgumentError("all coordinate arrays must have the same size"))
  end
end
