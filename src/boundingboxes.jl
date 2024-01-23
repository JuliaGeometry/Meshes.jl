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

function boundingbox(r::Ray{Dim,T}) where {Dim,T}
  lower(p, v) = v < 0 ? typemin(T) : p
  upper(p, v) = v > 0 ? typemax(T) : p
  p = r(0)
  v = r(1) - r(0)
  l = lower.(coordinates(p), v)
  u = upper.(coordinates(p), v)
  Box(Point(l), Point(u))
end

function boundingbox(s::Sphere{Dim,T}) where {Dim,T}
  c = center(s)
  r = radius(s)
  r⃗ = Vec(ntuple(i -> r, Dim))
  Box(c - r⃗, c + r⃗)
end

function boundingbox(p::ParaboloidSurface{T}) where {T}
  v = apex(p)
  r = radius(p)
  f = focallength(p)
  Box(v + Vec(-r, -r, T(0)), v + Vec(r, r, r^2 / (4f)))
end

boundingbox(t::Torus) = _pboxes(pointify(t))

boundingbox(g::CartesianGrid) = Box(extrema(g)...)

boundingbox(g::RectilinearGrid) = Box(extrema(g)...)

boundingbox(g::TransformedGrid{Dim,T,<:CartesianGrid{Dim,T}}) where {Dim,T} =
  boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(g::TransformedGrid{Dim,T,<:RectilinearGrid{Dim,T}}) where {Dim,T} =
  boundingbox(parent(g)) |> transform(g) |> boundingbox

boundingbox(m::Mesh) = _pboxes(vertices(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

_bboxes(boxes) = _pboxes(point for box in boxes for point in extrema(box))

function _pboxes(points)
  p = first(points)
  T = coordtype(p)
  Dim = embeddim(p)
  xmin = MVector(ntuple(i -> typemax(T), Dim))
  xmax = MVector(ntuple(i -> typemin(T), Dim))
  for p in points
    x = coordinates(p)
    @. xmin = min(x, xmin)
    @. xmax = max(x, xmax)
  end
  Box(Point(xmin), Point(xmax))
end
