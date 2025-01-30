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

    RegularGrid(dims)
    RegularGrid(dim1, dim2, ...)

Finally, a Cartesian grid can be constructed by only passing the dimensions
`dims` as a tuple, or by passing each dimension `dim1`, `dim2`, ... separately.
In this case, the origin and spacing default to (0,0,...) and (1,1,...).

## Examples

Create a 3D grid with 100x100x50 hexahedrons:

```julia
julia> RegularGrid(100, 100, 50)
```

Create a 2D grid with 100 x 100 quadrangles and origin at (10.0, 20.0):

```julia
julia> RegularGrid((100, 100), (10.0, 20.0), (1.0, 1.0))
```

Create a 1D grid from -1 to 1 with 100 segments:

```julia
julia> RegularGrid((-1.0,), (1.0,), dims=(100,))
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
  _checkorigin(origin)

  nc = CoordRefSystems.ncoords(C)

  if N ≠ nc
    throw(ArgumentError("""
    A $N-dimensional regular grid requires an origin with $N coordinates.
    The provided origin has $nc coordinates.
    """))
  end

  spac = _spacing(origin, spacing)

  RegularGrid{M,C,N,typeof(spac)}(origin, spac, offset, topology)
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

RegularGrid(
  dims::Dims{Dim},
  origin::NTuple{Dim,Number},
  spacing::NTuple{Dim,Number},
  offset::Dims{Dim}=ntuple(i -> 1, Dim)
) where {Dim} = RegularGrid(dims, Point(origin), spacing, offset)

function RegularGrid(start::Point, finish::Point, spacing::NTuple{N,Number}) where {N}
  _checkorigin(start)
  svals, fvals = _startfinish(start, finish)
  spac = _spacing(start, spacing)
  dims = ceil.(Int, (fvals .- svals) ./ spac)
  RegularGrid(dims, start, spac)
end

RegularGrid(start::NTuple{Dim,Number}, finish::NTuple{Dim,Number}, spacing::NTuple{Dim,Number}) where {Dim} =
  RegularGrid(Point(start), Point(finish), spacing)

function RegularGrid(start::Point, finish::Point; dims::Dims=ntuple(i -> 100, CoordRefSystems.ncoords(crs(start))))
  _checkorigin(start)
  svals, fvals = _startfinish(start, finish)
  spacing = (fvals .- svals) ./ dims
  RegularGrid(dims, start, spacing)
end

RegularGrid(start::NTuple{Dim,Number}, finish::NTuple{Dim,Number}; dims::Dims{Dim}=ntuple(i -> 100, Dim)) where {Dim} =
  RegularGrid(Point(start), Point(finish); dims)

function RegularGrid(dims::Dims{Dim}) where {Dim}
  origin = ntuple(i -> 0.0, Dim)
  spacing = ntuple(i -> 1.0, Dim)
  offset = ntuple(i -> 1, Dim)
  RegularGrid(dims, origin, spacing, offset)
end

RegularGrid(dims::Int...) = RegularGrid(dims)

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

function _checkorigin(origin)
  if manifold(origin) <: 🌐 && !(crs(origin) <: LatLon)
    throw(ArgumentError("regular spacing on `🌐` requires `LatLon` coordinates"))
  end
end

function _spacing(origin, spacing)
  C = crs(origin)
  T = CoordRefSystems.mactype(C)
  nc = CoordRefSystems.ncoords(C)
  us = CoordRefSystems.units(C)
  ntuple(i -> numconvert(T, withunit(spacing[i], us[i])), nc)
end

function _startfinish(start::Point{<:𝔼}, finish::Point{<:𝔼})
  scoords = coords(start)
  fcoords = convert(crs(start), coords(finish))
  svals = CoordRefSystems.values(scoords)
  fvals = CoordRefSystems.values(fcoords)
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
