# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StructuredGrid(X, Y, Z, ...)
    StructuredGrid{M,C}(X, Y, Z, ...)

A structured grid with vertices at sorted coordinates `X`, `Y`, `Z`, ...,
manifold `M` (default to `𝔼`) and CRS type `C` (default to `Cartesian`).

## Examples

Create a 2D structured grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> X = repeat(0.0:0.2:1.0, 1, 6)
julia> Y = repeat([0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
julia> StructuredGrid(X, Y)
```
"""
struct StructuredGrid{M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractArray}} <: Grid{M,C,N}
  XYZ::X
  topology::GridTopology{N}
  StructuredGrid{M,C,N,X}(XYZ, topology) where {M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractArray}} = new(XYZ, topology)
end

function StructuredGrid{M,C}(XYZ::NTuple{N,AbstractArray}, topology::GridTopology{N}) where {M<:Manifold,C<:CRS,N}
  if M <: 🌐 && !(C <: LatLon)
    throw(ArgumentError("rectilinear grid on `🌐` requires `LatLon` coordinates"))
  end

  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)

  if N ≠ nc
    throw(ArgumentError("""
    A $N-dimensional structured grid requires a CRS with $N coordinates.
    The provided CRS has $nc coordinates.
    """))
  end

  XYZ′ = ntuple(i -> numconvert.(T, _withunit.(XYZ[i], us[i])), nc)
  StructuredGrid{M,C,N,typeof(XYZ′)}(XYZ′, topology)
end

function StructuredGrid{M,C}(XYZ::NTuple{N,AbstractArray}) where {M<:Manifold,C<:CRS,N}
  if !allequal(size(X) for X in XYZ)
    throw(ArgumentError("all coordinate arrays must be the same size"))
  end

  nd = ndims(first(XYZ))

  if nd ≠ N
    throw(ArgumentError("""
    A $N-dimensional structured grid requires coordinate arrays with $N dimensions.
    The provided coordinate arrays have $nd dimensions.
    """))
  end

  topology = GridTopology(size(first(XYZ)) .- 1)
  StructuredGrid{M,C}(XYZ, topology)
end

StructuredGrid{M,C}(XYZ::AbstractArray...) where {M<:Manifold,C<:CRS} = StructuredGrid{M,C}(XYZ)

function StructuredGrid(XYZ::NTuple{N,AbstractArray}) where {N}
  try
    L = promote_type(ntuple(i -> _lentype(eltype(XYZ[i])), N)...)
    M = 𝔼{N}
    C = Cartesian{NoDatum,N,L}
    StructuredGrid{M,C}(XYZ)
  catch
    throw(ArgumentError("invalid units for cartesian coordinates"))
  end
end

StructuredGrid(XYZ::AbstractArray...) = StructuredGrid(XYZ)

function vertex(g::StructuredGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  Point(ctor(ntuple(d -> g.XYZ[d][ijk...], paramdim(g))...))
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
