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

boundingbox(p::Polytope) = _pboxes(vertices(p))

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
  Box(Point(coords(l)), Point(coords(u)))
end

function boundingbox(s::Sphere{Dim}) where {Dim}
  c = center(s)
  r = radius(s)
  r⃗ = Vec(ntuple(i -> r, Dim))
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
  vs = (1,)
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

boundingbox(g::CartesianGrid) = Box(extrema(g)...)

boundingbox(g::RectilinearGrid) = Box(extrema(g)...)

boundingbox(g::TransformedGrid{Dim,<:CartesianGrid{Dim}}) where {Dim} =
  boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(g::TransformedGrid{Dim,<:RectilinearGrid{Dim}}) where {Dim} =
  boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(m::Mesh) = _pboxes(vertices(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

_bboxes(boxes) = _pboxes(point for box in boxes for point in extrema(box))

function _pboxes(points)
  p = first(points)
  ℒ = lentype(p)
  Dim = embeddim(p)
  xmin = MVector(ntuple(i -> typemax(ℒ), Dim))
  xmax = MVector(ntuple(i -> typemin(ℒ), Dim))
  for p in points
    x = to(p)
    @. xmin = min(x, xmin)
    @. xmax = max(x, xmax)
  end
  Box(Point(coords(xmin)), Point(coords(xmax)))
end
