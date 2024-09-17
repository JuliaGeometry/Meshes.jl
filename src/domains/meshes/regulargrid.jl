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

    RegularGrid(start, finish, dims=dims)

Alternatively, construct a regular grid from a `start` point to a `finish`
with dimensions `dims`.

    RegularGrid(start, finish, spacing)

Alternatively, construct a regular grid from a `start` point to a `finish`
point using a given `spacing`.

## Examples

```
RegularGrid((10, 20), Point(LatLon(30.0°, 60.0°)), (1.0, 1.0)) # add coordinate units to spacing
RegularGrid((10, 20), Point(Polar(0.0cm, 0.0rad)), (10.0mm, 1.0rad)) # convert spacing units to coordinate units
RegularGrid((10, 20), Point(Marcator(0.0, 0.0)), (1.5, 1.5))
RegularGrid((10, 20, 30), Point(Cylindrical(0.0, 0.0, 0.0)), (3.0, 2.0, 1.0))
```

See also [`CartesianGrid`](@ref).
"""
struct RegularGrid{M<:Manifold,C<:CRS,N,S<:NTuple{N,Quantity}} <: Grid{M,C,N}
  origin::Point{M,C}
  spacing::S
  offset::Dims{N}
  topology::GridTopology{N}

  function RegularGrid{M,C,N,S}(origin, spacing, offset, topology) where {M<:Manifold,C<:CRS,N,S<:NTuple{N,Quantity}}
    if !all(s -> s > zero(s), spacing)
      throw(ArgumentError("spacing must be positive"))
    end
    new(origin, spacing, offset, topology)
  end
end

function RegularGrid(
  origin::Point{M,C},
  spacing::NTuple{N,Number},
  offset::Dims{N},
  topology::GridTopology{N}
) where {M<:Manifold,C<:CRS,N}
  if M <: 🌐 && !(C <: LatLon)
    throw(ArgumentError("regular spacing on `🌐` requires `LatLon` coordinates"))
  end

  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)

  if N ≠ nc
    throw(ArgumentError("""
    A $N-dimensional regular grid requires an origin with $N coordinates.
    The provided origin has $nc coordinates.
    """))
  end

  sp = ntuple(i -> numconvert(T, withunit(spacing[i], us[i])), nc)

  RegularGrid{M,C,N,typeof(sp)}(origin, sp, offset, topology)
end

function RegularGrid(
  dims::Dims{N},
  origin::Point,
  spacing::NTuple{N,Number},
  offset::Dims{N}=ntuple(i -> 1, N)
) where {N}
  if !all(>(0), dims)
    throw(ArgumentError("dimensions must be positive"))
  end
  RegularGrid(origin, spacing, offset, GridTopology(dims))
end

function RegularGrid(start::Point, finish::Point; dims::Dims=ntuple(i -> 100, CoordRefSystems.ncoords(crs(start))))
  svals, fvals = _startfinish(start, finish)
  spacing = (fvals .- svals) ./ dims
  RegularGrid(dims, start, spacing)
end

function RegularGrid(start::Point, finish::Point, spacing::NTuple{N,Number}) where {N}
  svals, fvals = _startfinish(start, finish)
  dims = ceil.(Int, (fvals .- svals) ./ spacing)
  RegularGrid(dims, start, spacing)
end

spacing(g::RegularGrid) = g.spacing

offset(g::RegularGrid) = g.offset

function vertex(g::RegularGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  orig = CoordRefSystems.values(coords(g.origin))
  vals = orig .+ (ijk .- g.offset) .* g.spacing
  Point(ctor(vals...))
end

@generated function xyz(g::RegularGrid{M,C,N}) where {M,C,N}
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

function _startfinish(start::Point{<:𝔼}, finish::Point{<:𝔼})
  finish′ = Point(convert(crs(start), coords(finish)))
  svals = CoordRefSystems.values(coords(start))
  fvals = CoordRefSystems.values(coords(finish′))
  svals, fvals
end

function _startfinish(start::Point{<:🌐}, finish::Point{<:🌐})
  slatlon = convert(LatLon, coords(start))
  flatlon = convert(LatLon, coords(finish))
  slon = flatlon.lon < slatlon.lon ? slatlon.lon - 360u"°" : slatlon.lon
  svals = (slatlon.lat, slon)
  fvals = (flatlon.lat, flatlon.lon)
  svals, fvals
end
