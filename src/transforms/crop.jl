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
  box = preprocess(t, d)
  n = view(d, box)
  n, nothing
end

function apply(t::Crop, g::CartesianGrid)
  box = preprocess(t, g)
  range = cartesianrange(g, box)
  g[range], nothing
end

function apply(t::Crop, g::RectilinearGrid)
  box = preprocess(t, g)
  range = cartesianrange(g, box)
  g[range], nothing
end

function apply(t::Crop, g::Grid{𝔼{2}})
  box = preprocess(t, g)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  sx, sy = size(g)
  xstart = findfirst(1:sx) do i
    p = centroid(g[i, 1])
    c = convert(Cartesian, coords(p))
    c.x ≥ min.x
  end
  xstop = findlast(1:sx) do i
    p = centroid(g[i, 1])
    c = convert(Cartesian, coords(p))
    c.x ≤ max.x
  end
  ystart = findfirst(1:sy) do i
    p = centroid(g[1, i])
    c = convert(Cartesian, coords(p))
    c.y ≥ min.y
  end
  ystop = findlast(1:sy) do i
    p = centroid(g[1, i])
    c = convert(Cartesian, coords(p))
    c.y ≤ max.y
  end
  g[xstart:xstop, ystart:ystop], nothing
end

function apply(t::Crop, g::Grid{𝔼{3}})
  box = preprocess(t, g)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  sx, sy, sz = size(g)
  xstart = findfirst(1:sx) do i
    p = centroid(g[i, 1, 1])
    c = convert(Cartesian, coords(p))
    c.x ≥ min.x
  end
  xstop = findlast(1:sx) do i
    p = centroid(g[i, 1, 1])
    c = convert(Cartesian, coords(p))
    c.x ≤ max.x
  end
  ystart = findfirst(1:sy) do i
    p = centroid(g[1, i, 1])
    c = convert(Cartesian, coords(p))
    c.y ≥ min.y
  end
  ystop = findlast(1:sy) do i
    p = centroid(g[1, i, 1])
    c = convert(Cartesian, coords(p))
    c.y ≤ max.y
  end
  zstart = findfirst(1:sz) do i
    p = centroid(g[1, 1, i])
    c = convert(Cartesian, coords(p))
    c.z ≥ min.z
  end
  zstop = findlast(1:sz) do i
    p = centroid(g[1, 1, i])
    c = convert(Cartesian, coords(p))
    c.z ≤ max.z
  end
  g[xstart:xstop, ystart:ystop, zstart:zstop], nothing
end

function apply(t::Crop, g::Grid{🌐})
  box = preprocess(t, g)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  slon, slat = size(g)
  lonstart = findfirst(1:slon) do i
    p = centroid(g[i, 1])
    c = convert(LatLon, coords(p))
    c.lon ≥ min.lon
  end
  lonstop = findlast(1:slon) do i
    p = centroid(g[i, 1])
    c = convert(LatLon, coords(p))
    c.lon ≤ max.lon
  end
  latstart = findfirst(1:slat) do i
    p = centroid(g[1, i])
    c = convert(LatLon, coords(p))
    c.lat ≥ min.lat
  end
  latstop = findlast(1:slat) do i
    p = centroid(g[1, i])
    c = convert(LatLon, coords(p))
    c.lat ≤ max.lat
  end
  g[lonstart:lonstop, latstart:latstop], nothing
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
