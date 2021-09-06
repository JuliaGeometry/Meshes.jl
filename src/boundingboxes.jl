# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -----------
# GEOMETRIES
# -----------

"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(g::Geometry) = boundingbox(vertices(g))

boundingbox(p::Primitive) = boundingbox(boundary(p))

boundingbox(b::Box) = b

function boundingbox(s::Sphere{Dim,T}) where {Dim,T}
  c = center(s)
  r = radius(s)
  r⃗ = Vec(ntuple(i->r, Dim))
  Box(c - r⃗, c + r⃗)
end

# --------
# DOMAINS
# --------

"""
    boundingbox(domain)

Axis-aligned bounding box of the `domain`.
"""
boundingbox(domain::Domain) = boundingbox(collect(domain))

boundingbox(grid::CartesianGrid) = Box(extrema(grid)...)

# -----
# DATA
# -----

boundingbox(data::Data) = boundingbox(domain(data))

# --------
# VECTORS
# --------

function boundingbox(points::AbstractVector{Point{Dim,T}}) where {Dim,T}
  xmin = MVector(ntuple(i->typemax(T), Dim))
  xmax = MVector(ntuple(i->typemin(T), Dim))
  for p in points
    x = coordinates(p)
    @. xmin = min(x, xmin)
    @. xmax = max(x, xmax)
  end
  Box(Point(xmin), Point(xmax))
end

boundingbox(boxes::AbstractVector{<:Box{Dim}}) where {Dim} =
  boundingbox([point for box in boxes for point in extrema(box)])

boundingbox(geometries::AbstractVector{<:Geometry}) =
  boundingbox(boundingbox.(geometries))
