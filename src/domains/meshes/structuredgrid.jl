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
struct StructuredGrid{M<:Manifold,C<:CRS,N,NT<:NTuple{N,AbstractArray}} <: Grid{M,C,N}
  coords::NT
  topology::GridTopology{N}
  StructuredGrid{M,C,N,NT}(coords, topology) where {M<:Manifold,C<:CRS,N,NT<:NTuple{N,AbstractArray}} =
    new(coords, topology)
end

function StructuredGrid{M,C}(coords::NTuple{N,AbstractArray}, topology::GridTopology{N}) where {M<:Manifold,C<:CRS,N}
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

  coords′ = ntuple(i -> numconvert.(T, withunit.(coords[i], us[i])), nc)

  StructuredGrid{M,C,N,typeof(coords′)}(coords′, topology)
end

function StructuredGrid{M,C}(coords::NTuple{N,AbstractArray}) where {M<:Manifold,C<:CRS,N}
  if !allequal(size(X) for X in coords)
    throw(ArgumentError("all coordinate arrays must be the same size"))
  end

  nd = ndims(first(coords))

  if nd ≠ N
    throw(ArgumentError("""
    A $N-dimensional structured grid requires coordinate arrays with $N dimensions.
    The provided coordinate arrays have $nd dimensions.
    """))
  end

  topology = GridTopology(size(first(coords)) .- 1)
  StructuredGrid{M,C}(coords, topology)
end

StructuredGrid{M,C}(coords::AbstractArray...) where {M<:Manifold,C<:CRS} = StructuredGrid{M,C}(coords)

function StructuredGrid(coords::NTuple{N,AbstractArray}) where {N}
  L = promote_type(ntuple(i -> aslentype(eltype(coords[i])), N)...)
  M = 𝔼{N}
  C = Cartesian{NoDatum,N,L}
  StructuredGrid{M,C}(coords)
end

StructuredGrid(coords::AbstractArray...) = StructuredGrid(coords)

function vertex(g::StructuredGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  Point(ctor(ntuple(d -> g.coords[d][ijk...], paramdim(g))...))
end

coordarrays(g::StructuredGrid) = g.coords

@generated function Base.getindex(g::StructuredGrid{M,C,N}, I::CartesianIndices) where {M,C,N}
  exprs = ntuple(i -> :(g.coords[$i][cinds]), N)

  quote
    @boundscheck _checkbounds(g, I)
    dims = size(I)
    cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
    coords = ($(exprs...),)
    StructuredGrid{M,C}(coords, GridTopology(dims))
  end
end

function Base.summary(io::IO, g::StructuredGrid)
  join(io, size(g), "×")
  print(io, " StructuredGrid")
end
