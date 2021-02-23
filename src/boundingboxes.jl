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
boundingbox(geom::Geometry) = boundingbox(vertices(geom))

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
boundingbox(domain::Domain) = boundingbox(PointSet(domain))

boundingbox(pset::PointSet) = boundingbox([pset[i] for i in 1:nelements(pset)])

boundingbox(grid::CartesianGrid) = Box(extrema(grid)...)

# ---------
# FALLBACK
# ---------

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
