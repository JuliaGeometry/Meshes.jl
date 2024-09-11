# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RectilinearGrid(x, y, z, ...)
    RectilinearGrid{M,C}(x, y, z, ...)

A rectilinear grid with vertices at sorted coordinates `x`, `y`, `z`, ...,
manifold `M` (default to `𝔼`) and CRS type `C` (default to `Cartesian`).

## Examples

Create a 2D rectilinear grid with regular spacing in `x` dimension
and irregular spacing in `y` dimension:

```julia
julia> x = 0.0:0.2:1.0
julia> y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
julia> RectilinearGrid(x, y)
```
"""
struct RectilinearGrid{M<:Manifold,C<:CRS,N,NT<:NTuple{N,AbstractVector}} <: Grid{M,C,N}
  coords::NT
  topology::GridTopology{N}
  RectilinearGrid{M,C,N,NT}(coords, topology) where {M<:Manifold,C<:CRS,N,NT<:NTuple{N,AbstractVector}} =
    new(coords, topology)
end

function RectilinearGrid{M,C}(coords::NTuple{N,AbstractVector}, topology::GridTopology{N}) where {M<:Manifold,C<:CRS,N}
  if M <: 🌐 && !(C <: LatLon)
    throw(ArgumentError("rectilinear grid on `🌐` requires `LatLon` coordinates"))
  end

  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)

  if N ≠ nc
    throw(ArgumentError("""
    A $N-dimensional rectilinear grid requires a CRS with $N coordinates.
    The provided CRS has $nc coordinates.
    """))
  end

  coords′ = ntuple(i -> numconvert.(T, withunit.(coords[i], us[i])), nc)

  RectilinearGrid{M,C,N,typeof(coords′)}(coords′, topology)
end

function RectilinearGrid{M,C}(coords::NTuple{N,AbstractVector}) where {M<:Manifold,C<:CRS,N}
  topology = GridTopology(length.(coords) .- 1)
  RectilinearGrid{M,C}(coords, topology)
end

RectilinearGrid{M,C}(coords::AbstractVector...) where {M<:Manifold,C<:CRS} = RectilinearGrid{M,C}(coords)

function RectilinearGrid(coords::NTuple{N,AbstractVector}) where {N}
  L = promote_type(ntuple(i -> aslentype(eltype(coords[i])), N)...)
  M = 𝔼{N}
  C = Cartesian{NoDatum,N,L}
  RectilinearGrid{M,C}(coords)
end

RectilinearGrid(coords::AbstractVector...) = RectilinearGrid(coords)

function vertex(g::RectilinearGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  Point(ctor(getindex.(g.coords, ijk)...))
end

coordvectors(g::RectilinearGrid) = g.coords

coordarrays(g::RectilinearGrid) = coordarrays(coordvectors(g))

@generated function Base.getindex(g::RectilinearGrid{M,C,N}, I::CartesianIndices) where {M,C,N}
  exprs = ntuple(N) do i
    :(g.coords[$i][start[$i]:stop[$i]])
  end

  quote
    @boundscheck _checkbounds(g, I)
    dims = size(I)
    start = Tuple(first(I))
    stop = Tuple(last(I)) .+ 1
    coords = ($(exprs...),)
    RectilinearGrid{M,C}(coords, GridTopology(dims))
  end
end

function Base.summary(io::IO, g::RectilinearGrid)
  join(io, size(g), "×")
  print(io, " RectilinearGrid")
end
