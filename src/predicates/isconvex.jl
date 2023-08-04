# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isconvex(geometry)

Tells whether or not the `geometry` is convex.
"""
function isconvex end

isconvex(::Point) = true

isconvex(::Segment) = true

isconvex(::Ray) = true

isconvex(::Line) = true

isconvex(::Plane) = true

isconvex(::Box) = true

isconvex(::Ball) = true

isconvex(::Sphere) = false

isconvex(::Disk) = true

isconvex(::Circle) = false

isconvex(::Cone) = true

isconvex(::ConeSurface) = false

isconvex(::Cylinder) = true

isconvex(::CylinderSurface) = false

isconvex(::Torus) = false

isconvex(::Triangle) = true

isconvex(::Tetrahedron) = true

isconvex(p::Polygon{2}) = Set(vertices(convexhull(p))) == Set(vertices(p))

isconvex(p::Polygon{3}) = isconvex(proj2D(p))

isconvex(m::Multi) = Set(vertices(convexhull(m))) == Set(vertices(m))

# --------------
# OPTIMIZATIONS
# --------------

function isconvex(q::Quadrangle{2})
  v = vertices(q)
  d1 = Segment(v[1], v[3])
  d2 = Segment(v[2], v[4])
  intersects(d1, d2)
end

function isconvex(q::Quadrangle{3})
  v = vertices(q)
  iscoplanar(v...) && isconvex(proj2D(q))
end
