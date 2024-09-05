# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Crop(x=(xmin, xmax), y=(ymin, ymax), z=(zmin, zmax))
    Crop(lat=(latmin, latmax), lon=(lonmin, lonmax))

Retain the domain geometries within `x` limits [`xmax`,`xmax`],
`y` limits [`ymax`,`ymax`] and `z` limits [`zmin`,`zmax`] in length units
(default to meters), or within latitude limits [`latmin`,`latmax`]
and longitude limits [`lonmin`,`lonmax`] in degree units.

## Examples

```julia
Crop(x=(2, 4))
Crop(x=(1u"km", 3u"km"))
Crop(y=(1.2, 1.8), z=(2.4, 3.0))
Crop(lat=(30, 60))
Crop(lon=(45u"°", 90u"°"))
```
"""
struct Crop{T} <: GeometricTransform
  limits::T
end

Crop(; kwargs...) = Crop(values(kwargs))

parameters(t::Crop) = (; limits=t.limits)

preprocess(t::Crop, d::Domain) = _crop(boundingbox(d), t.limits)

function apply(t::Crop, d::Domain)
  b = preprocess(t, d)
  view(d, b), nothing
end

function apply(t::Crop, g::CartesianGrid)
  b = preprocess(t, g)
  g[cartesianrange(g, b)], nothing
end

function apply(t::Crop, g::RectilinearGrid)
  b = preprocess(t, g)
  g[cartesianrange(g, b)], nothing
end

function apply(t::Crop, g::Grid{𝔼{2}})
  box = preprocess(t, g)
  bmin = convert(Cartesian, coords(minimum(box)))
  bmax = convert(Cartesian, coords(maximum(box)))
  nx, ny = vsize(g)

  gmin = convert(Cartesian, coords(vertex(g, (1, 1))))
  gxmin = gmin.x
  gymin = gmin.y
  gxmax = convert(Cartesian, coords(vertex(g, (nx, 1)))).x
  gymax = convert(Cartesian, coords(vertex(g, (1, ny)))).y

  xmin = max(bmin.x, gxmin)
  xmax = min(bmax.x, gxmax)
  ymin = max(bmin.y, gymin)
  ymax = min(bmax.y, gymax)

  iₛ = findlast(1:nx) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ xmin
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ xmax
  end
  jₛ = findlast(1:ny) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≤ ymin
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≥ ymax
  end
  
  if iₛ == iₑ || jₛ == jₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  g[iₛ:(iₑ - 1), jₛ:(jₑ - 1)], nothing
end

function apply(t::Crop, g::Grid{𝔼{3}})
  box = preprocess(t, g)
  bmin = convert(Cartesian, coords(minimum(box)))
  bmax = convert(Cartesian, coords(maximum(box)))
  nx, ny, nz = vsize(g)

  gmin = convert(Cartesian, coords(vertex(g, (1, 1, 1))))
  gxmin = gmin.x
  gymin = gmin.y
  gzmin = gmin.z
  gxmax = convert(Cartesian, coords(vertex(g, (nx, 1, 1)))).x
  gymax = convert(Cartesian, coords(vertex(g, (1, ny, 1)))).y
  gzmax = convert(Cartesian, coords(vertex(g, (1, 1, nz)))).z

  xmin = max(bmin.x, gxmin)
  xmax = min(bmax.x, gxmax)
  ymin = max(bmin.y, gymin)
  ymax = min(bmax.y, gymax)
  zmin = max(bmin.z, gzmin)
  zmax = min(bmax.z, gzmax)

  iₛ = findlast(1:nx) do i
    p = vertex(g, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ xmin
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(g, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ xmax
  end
  jₛ = findlast(1:ny) do i
    p = vertex(g, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≤ ymin
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(g, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≥ ymax
  end
  kₛ = findlast(1:nz) do i
    p = vertex(g, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≤ zmin
  end
  kₑ = findfirst(1:nz) do i
    p = vertex(g, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≥ zmax
  end

  if iₛ == iₑ || jₛ == jₑ || kₛ == kₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  g[iₛ:(iₑ - 1), jₛ:(jₑ - 1), kₛ:(kₑ - 1)], nothing
end

function apply(t::Crop, g::Grid{🌐})
  box = preprocess(t, g)
  bmin = convert(LatLon, coords(minimum(box)))
  bmax = convert(LatLon, coords(maximum(box)))
  nlon, nlat = vsize(g)

  gmin = convert(Cartesian, coords(vertex(g, (1, 1))))
  glonmin = gmin.lon
  glatmin = gmin.lat
  glonmax = convert(Cartesian, coords(vertex(g, (nlon, 1)))).lon
  glatmax = convert(Cartesian, coords(vertex(g, (1, nlat)))).lat

  lonmin = max(bmin.lon, glonmin)
  lonmax = min(bmax.lon, glonmax)
  latmin = max(bmin.lat, glatmin)
  latmax = min(bmax.lat, glatmax)

  iₛ = findlast(1:nlon) do i
    p = vertex(g, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≤ lonmin
  end
  iₑ = findfirst(1:nlon) do i
    p = vertex(g, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≥ lonmax
  end
  jₛ = findlast(1:nlat) do i
    p = vertex(g, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≤ latmin
  end
  jₑ = findfirst(1:nlat) do i
    p = vertex(g, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≥ latmax
  end

  if iₛ == iₑ || jₛ == jₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  g[iₛ:(iₑ - 1), jₛ:(jₑ - 1)], nothing
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _crop(box::Box{<:𝔼}, limits)
  lims = _xyzlimits(limits)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  xyzmin, xyzmax = _xyzminmax(min, max, lims)
  Box(withcrs(box, xyzmin), withcrs(box, xyzmax))
end

function _crop(box::Box{🌐}, limits)
  lims = _latlonlimits(limits)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  latmin, latmax = isnothing(lims.lat) ? (min.lat, max.lat) : lims.lat
  lonmin, lonmax = isnothing(lims.lon) ? (min.lon, max.lon) : lims.lon
  Box(withcrs(box, (latmin, lonmin), LatLon), withcrs(box, (latmax, lonmax), LatLon))
end

_xyzlimits(limits) = (
  x=haskey(limits, :x) ? _aslen.(limits.x) : nothing,
  y=haskey(limits, :y) ? _aslen.(limits.y) : nothing,
  z=haskey(limits, :z) ? _aslen.(limits.z) : nothing
)

_latlonlimits(limits) =
  (lat=haskey(limits, :lat) ? _asdeg.(limits.lat) : nothing, lon=haskey(limits, :lon) ? _asdeg.(limits.lon) : nothing)

function _xyzminmax(min::Cartesian2D, max::Cartesian2D, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  ymin, ymax = isnothing(lims.y) ? (min.y, max.y) : lims.y
  (xmin, ymin), (xmax, ymax)
end

function _xyzminmax(min::Cartesian3D, max::Cartesian3D, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  ymin, ymax = isnothing(lims.y) ? (min.y, max.y) : lims.y
  zmin, zmax = isnothing(lims.z) ? (min.z, max.z) : lims.z
  (xmin, ymin, zmin), (xmax, ymax, zmax)
end

_aslen(x::Len) = float(x)
_aslen(x::Number) = float(x) * u"m"
_aslen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

_asdeg(x::Deg) = float(x)
_asdeg(x::Number) = float(x) * u"°"
_asdeg(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))
