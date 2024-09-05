# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Crop(x=(xmin, xmax), y=(ymin, ymax), z=(zmin, zmax))
    Crop(lat=(latmin, latmax), lon=(lonmin, lonmax))

Retain the grid elements within `x` limits [`xmax`,`xmax`],
`y` limits [`ymax`,`ymax`] and `z` limits [`zmin`,`zmax`]
in length units (default to meters), or within `lat` limits
[`latmin`,`latmax`] and `lon` limits [`lonmin`,`lonmax`]
in degree units.

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

preprocess(t::Crop, g::Grid) = cartesianrange(g, _fixlimits(boundingbox(g), t.limits))

apply(t::Crop, g::Grid) = g[preprocess(t, g)], nothing

# -----------------
# HELPER FUNCTIONS
# -----------------

function _fixlimits(box::Box{<:𝔼}, limits)
  lims = _xyzlimits(limits)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  _minmax(min, max, lims)
end

function _fixlimits(box::Box{🌐}, limits)
  lims = _latlonlimits(limits)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  _minmax(min, max, lims)
end

_xyzlimits(limits) = (
  x=haskey(limits, :x) ? _aslen.(limits.x) : nothing,
  y=haskey(limits, :y) ? _aslen.(limits.y) : nothing,
  z=haskey(limits, :z) ? _aslen.(limits.z) : nothing
)

_latlonlimits(limits) =
  (lat=haskey(limits, :lat) ? _asdeg.(limits.lat) : nothing, lon=haskey(limits, :lon) ? _asdeg.(limits.lon) : nothing)

function _minmax(min::Cartesian2D, max::Cartesian2D, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  ymin, ymax = isnothing(lims.y) ? (min.y, max.y) : lims.y
  (xmin, xmax), (ymin, ymax)
end

function _minmax(min::Cartesian3D, max::Cartesian3D, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  ymin, ymax = isnothing(lims.y) ? (min.y, max.y) : lims.y
  zmin, zmax = isnothing(lims.z) ? (min.z, max.z) : lims.z
  (xmin, xmax), (ymin, ymax), (zmin, zmax)
end

function _minmax(min::LatLon, max::LatLon, lims)
  lonmin, lonmax = isnothing(lims.lon) ? (min.lon, max.lon) : lims.lon
  latmin, latmax = isnothing(lims.lat) ? (min.lat, max.lat) : lims.lat
  (lonmin, lonmax), (latmin, latmax)
end

_aslen(x::Len) = float(x)
_aslen(x::Number) = float(x) * u"m"
_aslen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

_asdeg(x::Deg) = float(x)
_asdeg(x::Number) = float(x) * u"°"
_asdeg(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))
