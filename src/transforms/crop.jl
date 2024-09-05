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
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  nx, ny = vsize(g)

  # check limits
  ivalid = any(1:nx) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    min.x ≤ c.x ≤ max.x
  end
  jvalid = any(1:ny) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    min.y ≤ c.y ≤ max.y
  end
  if !ivalid || !jvalid
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  iₛ = findlast(1:nx) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ min.x
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ max.x
  end
  jₛ = findlast(1:ny) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≤ min.y
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≥ max.y
  end
  irange = _fixindex(iₛ, 1):(_fixindex(iₑ, nx) - 1)
  jrange = _fixindex(jₛ, 1):(_fixindex(jₑ, ny) - 1)
  g[irange, jrange], nothing
end

function apply(t::Crop, g::Grid{𝔼{3}})
  box = preprocess(t, g)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  nx, ny, nz = vsize(g)

  # check limits
  ivalid = any(1:nx) do i
    p = vertex(g, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    min.x ≤ c.x ≤ max.x
  end
  jvalid = any(1:ny) do i
    p = vertex(g, (1, i, 1))
    c = convert(Cartesian, coords(p))
    min.y ≤ c.y ≤ max.y
  end
  kvalid = any(1:nz) do i
    p = vertex(g, (1, 1, i))
    c = convert(Cartesian, coords(p))
    min.z ≤ c.z ≤ max.z
  end
  if !ivalid || !jvalid || !kvalid
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  iₛ = findlast(1:nx) do i
    p = vertex(g, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ min.x
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(g, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ max.x
  end
  jₛ = findlast(1:ny) do i
    p = vertex(g, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≤ min.y
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(g, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≥ max.y
  end
  kₛ = findlast(1:nz) do i
    p = vertex(g, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≤ min.z
  end
  kₑ = findfirst(1:nz) do i
    p = vertex(g, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≥ max.z
  end
  irange = _fixindex(iₛ, 1):(_fixindex(iₑ, nx) - 1)
  jrange = _fixindex(jₛ, 1):(_fixindex(jₑ, ny) - 1)
  krange = _fixindex(kₛ, 1):(_fixindex(kₑ, nz) - 1)
  g[irange, jrange, krange], nothing
end

function apply(t::Crop, g::Grid{🌐})
  box = preprocess(t, g)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  nlon, nlat = vsize(g)

  # check limits
  ivalid = any(1:nlon) do i
    p = vertex(g, (i, 1))
    c = convert(Cartesian, coords(p))
    min.lon ≤ c.lon ≤ max.lon
  end
  jvalid = any(1:nlat) do i
    p = vertex(g, (1, i))
    c = convert(Cartesian, coords(p))
    min.lat ≤ c.lat ≤ max.lat
  end
  if !ivalid || !jvalid
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  iₛ = findlast(1:nlon) do i
    p = vertex(g, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≤ min.lon
  end
  iₑ = findfirst(1:nlon) do i
    p = vertex(g, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≥ max.lon
  end
  jₛ = findlast(1:nlat) do i
    p = vertex(g, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≤ min.lat
  end
  jₑ = findfirst(1:nlat) do i
    p = vertex(g, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≥ max.lat
  end
  irange = _fixindex(iₛ, 1):(_fixindex(iₑ, nlon) - 1)
  jrange = _fixindex(jₛ, 1):(_fixindex(jₑ, nlat) - 1)
  g[irange, jrange], nothing
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

_fixindex(i, default) = isnothing(i) ? default : i

_aslen(x::Len) = float(x)
_aslen(x::Number) = float(x) * u"m"
_aslen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

_asdeg(x::Deg) = float(x)
_asdeg(x::Number) = float(x) * u"°"
_asdeg(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))
