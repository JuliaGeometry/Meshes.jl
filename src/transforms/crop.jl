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
struct Crop{X,Y,Z} <: GeometricTransform
  x::X
  y::Y
  z::Z
end

Crop(; x=nothing, y=nothing, z=nothing) =
  Crop(isnothing(x) ? x : _aslen.(x), isnothing(y) ? y : _aslen.(y), isnothing(z) ? z : _aslen.(z))

parameters(t::Crop) = (; x=t.x, y=t.y, z=t.z)

function preprocess(t::Crop, d::Domain)
  bbox = boundingbox(d)
  bbox₁ = _overlaps(1, t.x, bbox)
  bbox₂ = _overlaps(2, t.y, bbox₁)
  bbox₃ = _overlaps(3, t.z, bbox₂)
  bbox₃
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
