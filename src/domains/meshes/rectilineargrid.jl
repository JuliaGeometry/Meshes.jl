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
struct RectilinearGrid{M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractVector}} <: Grid{M,C,N}
  xyz::X
  topology::GridTopology{N}
  RectilinearGrid{M,C,N,X}(xyz, topology) where {M<:Manifold,C<:CRS,N,X<:NTuple{N,AbstractVector}} = new(xyz, topology)
end

function RectilinearGrid{M,C}(xyz::NTuple{N,AbstractVector}, topology::GridTopology{N}) where {M<:Manifold,C<:CRS,N}
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

  xyz′ = ntuple(i -> numconvert.(T, withunit.(xyz[i], us[i])), nc)
  RectilinearGrid{M,C,N,typeof(xyz′)}(xyz′, topology)
end

function RectilinearGrid{M,C}(xyz::NTuple{N,AbstractVector}) where {M<:Manifold,C<:CRS,N}
  topology = GridTopology(length.(xyz) .- 1)
  RectilinearGrid{M,C}(xyz, topology)
end

RectilinearGrid{M,C}(xyz::AbstractVector...) where {M<:Manifold,C<:CRS} = RectilinearGrid{M,C}(xyz)

function RectilinearGrid(xyz::NTuple{N,AbstractVector}) where {N}
  L = promote_type(ntuple(i -> aslentype(eltype(xyz[i])), N)...)
  M = 𝔼{N}
  C = Cartesian{NoDatum,N,L}
  RectilinearGrid{M,C}(xyz)
end

RectilinearGrid(xyz::AbstractVector...) = RectilinearGrid(xyz)

function vertex(g::RectilinearGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  Point(ctor(getindex.(g.xyz, ijk)...))
end

xyz(g::RectilinearGrid) = g.xyz

XYZ(g::RectilinearGrid) = XYZ(xyz(g))

@generated function Base.getindex(g::RectilinearGrid{M,C,N}, I::CartesianIndices) where {M,C,N}
  exprs = ntuple(N) do i
    :(g.xyz[$i][start[$i]:stop[$i]])
  end

  quote
    @boundscheck _checkbounds(g, I)
    dims = size(I)
    start = Tuple(first(I))
    stop = Tuple(last(I)) .+ 1
    xyz = ($(exprs...),)
    RectilinearGrid{M,C}(xyz, GridTopology(dims))
  end
end

function Base.summary(io::IO, g::RectilinearGrid)
  join(io, size(g), "×")
  print(io, " RectilinearGrid")
end
