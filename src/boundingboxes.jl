# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    boundingbox(object)

Axis-aligned bounding box of `object`.
"""
function boundingbox end

# ----------
# FALLBACKS
# ----------

boundingbox(p::Polytope) = _pboxes(eachvertex(p))

boundingbox(p::Primitive) = boundingbox(boundary(p))

boundingbox(m::Multi) = _bboxes(boundingbox(g) for g in parent(m))

boundingbox(geoms) = _bboxes(boundingbox(g) for g in geoms)

# ----------------
# SPECIALIZATIONS
# ----------------

boundingbox(p::Point) = Box(p, p)

boundingbox(b::Box) = b

function boundingbox(r::Ray)
  lower(p, v) = v < zero(v) ? typemin(p) : p
  upper(p, v) = v > zero(v) ? typemax(p) : p
  p = r(0)
  v = r(1) - r(0)
  l = lower.(to(p), v)
  u = upper.(to(p), v)
  Box(withcrs(r, l), withcrs(r, u))
end

function boundingbox(s::Sphere)
  c = center(s)
  r = radius(s)
  r⃗ = Vec(ntuple(i -> r, embeddim(s)))
  Box(c - r⃗, c + r⃗)
end

function boundingbox(c::CylinderSurface)
  us = (0, 1 / 4, 1 / 2, 3 / 4)
  vs = (0, 1 / 2, 1)
  ps = [c(u, v) for (u, v) in Iterators.product(us, vs)]
  boundingbox(ps)
end

function boundingbox(c::ConeSurface)
  us = (0, 1 / 4, 1 / 2, 3 / 4)
  vs = (0,)
  ps = [c(u, v) for (u, v) in Iterators.product(us, vs)]
  boundingbox([ps; apex(c)])
end

function boundingbox(p::ParaboloidSurface)
  v = apex(p)
  r = radius(p)
  f = focallength(p)
  Box(v + Vec(-r, -r, zero(r)), v + Vec(r, r, r^2 / (4f)))
end

boundingbox(t::Torus) = _pboxes(pointify(t))

boundingbox(g::OrthoRegularGrid) = Box(extrema(g)...)

boundingbox(g::OrthoRectilinearGrid) = Box(extrema(g)...)

boundingbox(g::TransformedGrid{<:Any,<:Any,<:OrthoRegularGrid}) = boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(g::TransformedGrid{<:Any,<:Any,<:OrthoRectilinearGrid}) =
  boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(m::Mesh) = _pboxes(eachvertex(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

_bboxes(boxes) = _pboxes(point for box in boxes for point in extrema(box))

_pboxes(points) = _pboxes(manifold(first(points)), points)

function _pboxes(::Type{𝔼{1}}, points)
  p = first(points)
  ℒ = lentype(p)
  xmin = typemax(ℒ)
  xmax = typemin(ℒ)
  for p in points
    c = convert(Cartesian, coords(p))
    xmin = min(c.x, xmin)
    xmax = max(c.x, xmax)
  end
  Box(withcrs(p, (xmin,)), withcrs(p, (xmax,)))
end

function _pboxes(::Type{𝔼{2}}, points)
  p = first(points)
  ℒ = lentype(p)
  xmin, ymin = typemax(ℒ), typemax(ℒ)
  xmax, ymax = typemin(ℒ), typemin(ℒ)
  for p in points
    c = convert(Cartesian, coords(p))
    xmin = min(c.x, xmin)
    ymin = min(c.y, ymin)
    xmax = max(c.x, xmax)
    ymax = max(c.y, ymax)
  end
  Box(withcrs(p, (xmin, ymin)), withcrs(p, (xmax, ymax)))
end

function _pboxes(::Type{𝔼{3}}, points)
  p = first(points)
  ℒ = lentype(p)
  xmin, ymin, zmin = typemax(ℒ), typemax(ℒ), typemax(ℒ)
  xmax, ymax, zmax = typemin(ℒ), typemin(ℒ), typemin(ℒ)
  for p in points
    c = convert(Cartesian, coords(p))
    xmin = min(c.x, xmin)
    ymin = min(c.y, ymin)
    zmin = min(c.z, zmin)
    xmax = max(c.x, xmax)
    ymax = max(c.y, ymax)
    zmax = max(c.z, zmax)
  end
  Box(withcrs(p, (xmin, ymin, zmin)), withcrs(p, (xmax, ymax, zmax)))
end

function _pboxes(::Type{🌐}, points)
  p = first(points)
  T = numtype(lentype(p))
  lonmin, latmin = T(180) * u"°", T(90) * u"°"
  lonmax, latmax = T(-180) * u"°", T(-90) * u"°"
  for p in points
    c = convert(LatLon, coords(p))
    lonmin = min(c.lon, lonmin)
    latmin = min(c.lat, latmin)
    lonmax = max(c.lon, lonmax)
    latmax = max(c.lat, latmax)
  end
  Box(withcrs(p, (latmin, lonmin), LatLon), withcrs(p, (latmax, lonmax), LatLon))
end
