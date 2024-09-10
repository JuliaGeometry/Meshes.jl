# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularGrid(dims, origin, spacing)

A regular grid with dimensions `dims`, lower left corner at `origin`
and cell spacing `spacing`. The three arguments must have the same length.

    RegularGrid(dims, origin, spacing, offset)

A regular grid with dimensions `dims`, with lower left corner of element
`offset` at `origin` and cell spacing `spacing`.

## Examples

```
RegularGrid((10, 20), Point(LatLon(30.0°, 60.0°)), (1.0, 1.0)) # add coordinate units to spacing
RegularGrid((10, 20), Point(Polar(0.0cm, 0.0rad)), (10.0mm, 1.0rad)) # convert spacing units to coordinate units
RegularGrid((10, 20), Point(Marcator(0.0, 0.0)), (1.5, 1.5))
RegularGrid((10, 20, 30), Point(Cylindrical(0.0, 0.0, 0.0)), (3.0, 2.0, 1.0))
```

See also [`CartesianGrid`](@ref).
"""
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
  if manifold(origin) <: 🌐 && !(crs(origin) <: LatLon)
    throw(ArgumentError("regular spacing on `🌐` requires `LatLon` coordinates"))
  end

  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)
  ns = length(spacing)

  if N ≠ nc
    throw(ArgumentError("""
    A $N-dimensional regular grid requires an origin with $N coordinates.
    The provided origin has $nc coordinates.
    """))
  end

  if ns ≠ nc
    throw(ArgumentError("""
    A $N-dimensional regular grid requires $N spacing values.
    The provided spacing has $ns values.
    """))
  end

  sp = ntuple(i -> numconvert(T, _withunit(spacing[i], us[i])), nc)

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

@generated function xyz(g::RegularGrid{M,C,S,N}) where {M,C,S,N}
  exprs = ntuple(N) do i
    :(range(start=orig[$i], step=spac[$i], length=(dims[$i] + 1)))
  end

  quote
    dims = size(g)
    spac = spacing(g)
    orig = CoordRefSystems.values(coords(g.origin))
    ($(exprs...),)
  end
end

XYZ(g::RegularGrid) = XYZ(xyz(g))

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

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, g::RegularGrid)
  dims = join(size(g.topology), "×")
  name = prettyname(g)
  print(io, "$dims $name")
end

function Base.show(io::IO, ::MIME"text/plain", g::RegularGrid)
  summary(io, g)
  println(io)
  println(io, "├─ minimum: ", minimum(g))
  println(io, "├─ maximum: ", maximum(g))
  print(io, "└─ spacing: ", spacing(g))
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_withunit(x::Number, u) = x * u
_withunit(x::Quantity, u) = uconvert(u, x)
