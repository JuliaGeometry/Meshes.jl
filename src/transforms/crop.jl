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

function Crop(; x=nothing, y=nothing, z=nothing, lat=nothing, lon=nothing)
  if !isnothing(lat) || !isnothing(lon)
    Crop((lat=isnothing(lat) ? lat : _asdeg.(lat), lon=isnothing(lon) ? lon : _asdeg.(lon)))
  else
    Crop((x=isnothing(x) ? x : _aslen.(x), y=isnothing(y) ? y : _aslen.(y), z=isnothing(z) ? z : _aslen.(z)))
  end
end

parameters(t::Crop) = (; limits=t.limits)

function preprocess(t::Crop, d::Domain{<:𝔼})
  bbox = boundingbox(d)
  bbox₁ = _overlaps(1, t.limits.x, bbox)
  bbox₂ = _overlaps(2, t.limits.y, bbox₁)
  bbox₃ = _overlaps(3, t.limits.z, bbox₂)
  bbox₃
end

function preprocess(t::Crop, d::Domain{<:🌐})
  bbox = boundingbox(d)
  bbox₁ = _overlapslat(t.limits.lat, bbox)
  bbox₂ = _overlapslon(t.limits.lon, bbox₁)
  bbox₂
end

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

_aslen(x::Len) = float(x)
_aslen(x::Number) = float(x) * u"m"
_aslen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

_asdeg(x::Deg) = float(x)
_asdeg(x::Number) = float(x) * u"°"
_asdeg(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

function _overlaps(dim, lims, bbox)
  Dim = embeddim(bbox)
  if Dim < dim || isnothing(lims)
    bbox
  else
    lmin, lmax = lims
    min = to(minimum(bbox))
    max = to(maximum(bbox))
    nmin = Vec(ntuple(i -> i == dim ? lmin : min[i], Dim))
    nmax = Vec(ntuple(i -> i == dim ? lmax : max[i], Dim))
    Box(withcrs(bbox, nmin), withcrs(bbox, nmax))
  end
end

function _overlapslat(lims, bbox)
  point(lat, lon) = Point{🌐}(convert(crs(bbox), LatLon{datum(crs(bbox))}(lat, lon)))
  if isnothing(lims)
    bbox
  else
    lmin, lmax = lims
    min = convert(LatLon, coords(minimum(bbox)))
    max = convert(LatLon, coords(maximum(bbox)))
    nmin = point(lmin, min.lon)
    nmax = point(lmax, max.lon)
    Box(nmin, nmax)
  end
end

function _overlapslon(lims, bbox)
  point(lat, lon) = Point{🌐}(convert(crs(bbox), LatLon{datum(crs(bbox))}(lat, lon)))
  if isnothing(lims)
    bbox
  else
    lmin, lmax = lims
    min = convert(LatLon, coords(minimum(bbox)))
    max = convert(LatLon, coords(maximum(bbox)))
    nmin = point(min.lat, lmin)
    nmax = point(max.lat, lmax)
    Box(nmin, nmax)
  end
end
