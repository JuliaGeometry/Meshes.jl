# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Slice(x=(xmin, xmax), y=(ymin, ymax), z=(zmin, zmax))
    Slice(lat=(latmin, latmax), lon=(lonmin, lonmax))

Retain the domain elements within `x` limits [`xmax`,`xmax`],
`y` limits [`ymax`,`ymax`] and `z` limits [`zmin`,`zmax`]
in length units (default to meters), or within `lat` limits
[`latmin`,`latmax`] and `lon` limits [`lonmin`,`lonmax`]
in degree units.

## Examples

```julia
Slice(x=(1000km, 3000km))
Slice(x=(1000km, 2000km), y=(2000km, 5000km))
Slice(lon=(0°, 90°))
Slice(lon=(0°, 45°), lat=(0°, 45°))
```
"""
struct Slice{T} <: GeometricTransform
  limits::T
end

Slice(; kwargs...) = Slice(values(kwargs))

parameters(t::Slice) = (; limits=t.limits)

preprocess(t::Slice, d::Domain) = _sliceinds(d, _slicebox(boundingbox(d), t.limits))

apply(t::Slice, d::Domain) = _slice(d, preprocess(t, d)), nothing

# -----------------
# HELPER FUNCTIONS
# -----------------

_slice(d::Domain, inds) = view(d, inds)
_slice(g::Grid, inds::CartesianIndices) = getindex(g, inds)

_sliceinds(d::Domain, b) = indices(d, b)
_sliceinds(g::OrthoRegularGrid, b) = cartesianrange(g, b)
_sliceinds(g::OrthoRectilinearGrid, b) = cartesianrange(g, b)
_sliceinds(g::Grid{🌐}, b::Box{🌐}) = cartesianrange(g, b)

function _slicebox(box::Box{𝔼{2}}, limits)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  xmin, xmax = get(limits, :x, (min.x, max.x))
  ymin, ymax = get(limits, :y, (min.y, max.y))
  bmin = _aslen.((xmin, ymin))
  bmax = _aslen.((xmax, ymax))
  Box(withcrs(box, bmin), withcrs(box, bmax))
end

function _slicebox(box::Box{𝔼{3}}, limits)
  min = convert(Cartesian, coords(minimum(box)))
  max = convert(Cartesian, coords(maximum(box)))
  xmin, xmax = get(limits, :x, (min.x, max.x))
  ymin, ymax = get(limits, :y, (min.y, max.y))
  zmin, zmax = get(limits, :z, (min.z, max.z))
  bmin = _aslen.((xmin, ymin, zmin))
  bmax = _aslen.((xmax, ymax, zmax))
  Box(withcrs(box, bmin), withcrs(box, bmax))
end

function _slicebox(box::Box{🌐}, limits)
  min = convert(LatLon, coords(minimum(box)))
  max = convert(LatLon, coords(maximum(box)))
  latmin, latmax = get(limits, :lat, (min.lat, max.lat))
  lonmin, lonmax = get(limits, :lon, (min.lon, max.lon))
  bmin = _asdeg.((latmin, lonmin))
  bmax = _asdeg.((latmax, lonmax))
  Box(withcrs(box, bmin, LatLon), withcrs(box, bmax, LatLon))
end

_aslen(x::Len) = float(x)
_aslen(x::Number) = float(x) * u"m"
_aslen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

_asdeg(x::Deg) = float(x)
_asdeg(x::Number) = float(x) * u"°"
_asdeg(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))
