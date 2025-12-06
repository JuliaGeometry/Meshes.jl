# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularGrid(min, max; dims=dims)

A regular grid from `min` point to `max` point with dimensions `dims`.
The number of dimensions must match the number of coordinates of the points.

    RegularGrid(min, max, spacing)

Alternatively, construct a regular grid from `min` point to `max` point
by specifying the `spacing` for each dimension.

    RegularGrid(dims)
    RegularGrid(dim₁, dim₂, ...)

Alternatively, construct a regular grid with dimensions `dims = (dim₁, dim₂, ...)`,
min point at `(0m, 0m, ...)` and spacing equal to `(1m, 1m, ...)`.

    RegularGrid(origin, spacing, topology)

Finally, construct a regular grid with `origin` point, `spacing` and grid `topology`.
This method is available for advanced use cases involving periodic dimensions. See
[`GridTopology`](@ref) for more details.

## Examples

```julia
julia> RegularGrid((-1.0,), (1.0,), dims=(100,)) # 1D grid with 100 segments
julia> RegularGrid((0.0, 0.0), (10.0, 20.0), (1.0, 2.0)) # 2D grid with tall quadrangles
julia> RegularGrid(100, 100, 50) # 3D grid with 100x100x50 hexahedrons
```

See also [`CartesianGrid`](@ref).
"""
struct RegularGrid{M<:Manifold,C<:CRS,N,S<:NTuple{N,Quantity}} <: Grid{M,C,N}
  origin::Point{M,C}
  spacing::S
  topology::GridTopology{N}

  function RegularGrid{M,C,N,S}(origin, spacing, topology) where {M<:Manifold,C<:CRS,N,S<:NTuple{N,Quantity}}
    if !all(s -> s > zero(s), spacing)
      throw(ArgumentError("spacing must be positive"))
    end
    new(origin, spacing, topology)
  end
end

function RegularGrid(
  origin::Point{M,C},
  spacing::NTuple{N,Number},
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

  RegularGrid{M,C,N,typeof(spac)}(origin, spac, topology)
end

RegularGrid(origin::NTuple{N,Number}, spacing::NTuple{N,Number}, topology::GridTopology{N}) where {N} =
  RegularGrid(Point(origin), spacing, topology)

function RegularGrid(min::Point, max::Point, spacing::NTuple{N,Number}) where {N}
  _checkorigin(min)
  cmin, cmax = _minmaxcoords(min, max)
  spac = _spacing(min, spacing)
  dims = ceil.(Int, (cmax .- cmin) ./ spac)
  RegularGrid(min, spac, GridTopology(dims))
end

RegularGrid(min::NTuple{N,Number}, max::NTuple{N,Number}, spacing::NTuple{N,Number}) where {N} =
  RegularGrid(Point(min), Point(max), spacing)

function RegularGrid(min::Point, max::Point; dims::Dims=ntuple(i -> 100, CoordRefSystems.ncoords(crs(min))))
  _checkorigin(min)
  cmin, cmax = _minmaxcoords(min, max)
  spac = (cmax .- cmin) ./ dims
  RegularGrid(min, spac, GridTopology(dims))
end

RegularGrid(min::NTuple{N,Number}, max::NTuple{N,Number}; dims::Dims{N}=ntuple(i -> 100, N)) where {N} =
  RegularGrid(Point(min), Point(max); dims)

function RegularGrid(dims::Dims{N}) where {N}
  orig = ntuple(i -> 0.0, N)
  spac = ntuple(i -> 1.0, N)
  RegularGrid(orig, spac, GridTopology(dims))
end

RegularGrid(dims::Int...) = RegularGrid(dims)

spacing(g::RegularGrid) = g.spacing

function vertex(g::RegularGrid, ijk::Dims)
  ctor = CoordRefSystems.constructor(crs(g))
  orig = CoordRefSystems.values(coords(g.origin))
  vals = orig .+ (ijk .- 1) .* g.spacing
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
  orig = vertex(g, Tuple(first(I)))
  spac = spacing(g)
  topo = GridTopology(size(I), isperiodic(topology(g)))
  RegularGrid(orig, spac, topo)
end

function ==(g₁::RegularGrid, g₂::RegularGrid)
  orig₁ = CoordRefSystems.values(coords(g₁.origin))
  orig₂ = CoordRefSystems.values(coords(g₂.origin))
  orig₁ == orig₂ && g₁.spacing == g₂.spacing && g₁.topology == g₂.topology
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

function _minmaxcoords(min::Point{<:𝔼}, max::Point{<:𝔼})
  mincoords = coords(min)
  maxcoords = convert(crs(min), coords(max))
  minvalues = CoordRefSystems.values(mincoords)
  maxvalues = CoordRefSystems.values(maxcoords)
  minvalues, maxvalues
end

function _minmaxcoords(min::Point{<:🌐}, max::Point{<:🌐})
  mincoords = convert(LatLon, coords(min))
  maxcoords = convert(LatLon, coords(max))
  minvalues = (mincoords.lat, maxcoords.lon < mincoords.lon ? mincoords.lon - 360u"°" : mincoords.lon)
  maxvalues = (maxcoords.lat, maxcoords.lon)
  minvalues, maxvalues
end
