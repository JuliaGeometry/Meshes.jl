# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct RegularGrid{M<:Manifold,C<:CRS,S<:Tuple,N} <: Grid{M,C,N}
  origin::Point{M,C}
  spacing::S
  offset::Dims{N}
  topology::GridTopology{N}

  function RegularGrid{M,C,S,N}(origin, spacing, offset, topology) where {M<:Manifold,C<:CRS,S<:Tuple,N}
    if !all(s -> s > zero(s), spacing)
      throw(ArgumentError("spacing must be positive"))
    end
    new(origin, spacing, offset, topology)
  end
end

function RegularGrid(
  origin::Point{M,C},
  spacing::Tuple,
  offset::Dims{N},
  topology::GridTopology{N}
) where {M<:Manifold,C<:CRS,N}
  nc = CoordRefSystems.ncoords(C)

  if N ≠ nc
    throw(ArgumentError("the number of dimensions must be equal to the number of origin coordinates"))
  end

  if length(spacing) ≠ nc
    throw(ArgumentError("the number of spacings must be equal to the number of origin coordinates"))
  end

  us = CoordRefSystems.units(C)
  sp = ntuple(i -> float(_withunit(spacing[i], us[i])), nc)

  RegularGrid{M,C,typeof(sp),N}(origin, sp, offset, topology)
end

function RegularGrid(dims::Dims{N}, origin::Point, spacing::Tuple, offset::Dims{N}=ntuple(i -> 1, N)) where {N}
  if !all(>(0), dims)
    throw(ArgumentError("dimensions must be positive"))
  end
  RegularGrid(origin, spacing, offset, GridTopology(dims))
end

spacing(g::RegularGrid) = g.spacing

offset(g::RegularGrid) = g.offset

function vertex(g::RegularGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  orig = CoordRefSystems.values(coords(g.origin))
  vals = orig .+ (ijk .- g.offset) .* g.spacing
  Point(ctor(vals...))
end

function Base.getindex(g::RegularGrid, I::CartesianIndices)
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  offset = g.offset .- Tuple(first(I)) .+ 1
  RegularGrid(dims, g.origin, g.spacing, offset)
end

function ==(g₁::RegularGrid, g₂::RegularGrid)
  orig₁ = CoordRefSystems.values(coords(g₁.origin))
  orig₂ = CoordRefSystems.values(coords(g₂.origin))
  g₁.topology == g₂.topology && g₁.spacing == g₂.spacing && orig₁ .- orig₂ == (g₁.offset .- g₂.offset) .* g₁.spacing
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_withunit(x::Number, u) = x * u
_withunit(x::Quantity, u) = uconvert(u, x)
