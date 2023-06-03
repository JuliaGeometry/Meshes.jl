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

boundingbox(p::Polytope) = boundingbox(vertices(p))

boundingbox(p::Primitive) = boundingbox(boundary(p))

boundingbox(m::Multi) = boundingbox(collect(m))

boundingbox(d::Domain) = boundingbox(collect(d))

boundingbox(d::Data) = boundingbox(domain(d))

# ----------------
# SPECIALIZATIONS
# ----------------

boundingbox(p::Point) = Box(p, p)

boundingbox(b::Box) = b

function boundingbox(s::Sphere{Dim,T}) where {Dim,T}
  c = center(s)
  r = radius(s)
  r⃗ = Vec(ntuple(i -> r, Dim))
  Box(c - r⃗, c + r⃗)
end

boundingbox(g::Grid) = Box(extrema(g)...)

# ---------------
# IMPLEMENTATION
# ---------------

boundingbox(geoms::AbstractVector{<:Geometry}) = boundingbox(boundingbox.(geoms))

boundingbox(boxes::AbstractVector{<:Box{Dim}}) where {Dim} =
  boundingbox([point for box in boxes for point in extrema(box)])

function boundingbox(points::AbstractVector{Point{Dim,T}}) where {Dim,T}
  xmin = MVector(ntuple(i -> typemax(T), Dim))
  xmax = MVector(ntuple(i -> typemin(T), Dim))
  for p in points
    x = coordinates(p)
    @. xmin = min(x, xmin)
    @. xmax = max(x, xmax)
  end
  Box(Point(xmin), Point(xmax))
end
