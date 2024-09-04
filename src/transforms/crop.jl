# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Crop(x=(xmin, xmax), y=(ymin, ymax), z=(zmin, zmax))

Retain the domain geometries that intersect with `x` limits [`xmax`,`xmax`],
`y` limits [`ymax`,`ymax`] and `z` limits [`zmin`,`zmax`] in length units
(default to meters).

## Examples

```julia
Crop(x=(2, 4))
Crop(x=(1u"km", 3u"km"))
Crop(y=(1.2, 1.8), z=(2.4, 3.0))
```
"""
struct Crop{T} <: GeometricTransform
  limits::T
end

Crop(; kwargs...) = Crop(values(kwargs))

parameters(t::Crop) = (; limits=t.limits)

preprocess(t::Crop, d::Domain) = _cropbox(boundingbox(d), t.limits)

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

# -----------------
# HELPER FUNCTIONS
# -----------------

function _cropbox(box::Box{<:𝔼}, limits)
  lims = _xyzlimits(limits)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  xyzmin, xyzmax = _xyzminmax(manifold(box), min, max, lims)
  Box(withcrs(box, xyzmin), withcrs(box, xyzmax))
end

function _cropbox(box::Box{🌐}, limits)
  lims = _latlonlimits(limits)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  latmin, latmax = isnothing(lims.lat) ? (min.lat, max.lat) : lims.lat
  lonmin, lonmax = isnothing(lims.lon) ? (min.lon, max.lon) : lims.lon
  Box(withcrs(box, (latmin, lonmin), crs=LatLon), withcrs(box, (latmax, lonmax), crs=LatLon))
end

_xyzlimits(limits) = (
  x=haskey(limits, :x) ? _aslen.(limits.x) : nothing,
  y=haskey(limits, :y) ? _aslen.(limits.y) : nothing,
  z=haskey(limits, :z) ? _aslen.(limits.z) : nothing
)

_latlonlimits(limits) =
  (lat=haskey(limits, :lat) ? _asdeg.(limits.lat) : nothing, lon=haskey(limits, :lon) ? _asdeg.(limits.lon) : nothing)

function _xyzminmax(::Type{𝔼{1}}, min, max, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  (xmin,), (xmax,)
end

function _xyzminmax(::Type{𝔼{2}}, min, max, lims)
  xmin, xmax = isnothing(lims.x) ? (min.x, max.x) : lims.x
  ymin, ymax = isnothing(lims.y) ? (min.y, max.y) : lims.y
  (xmin, ymin), (xmax, ymax)
end

function _xyzminmax(::Type{𝔼{3}}, min, max, lims)
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
