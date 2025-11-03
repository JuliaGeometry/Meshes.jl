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

function isconvex(b::BezierCurve)
  if ncontrols(b) ≤ 2
    return true
  else
    ps = controls(b)
    p₁, p₂ = ps[begin], ps[begin + 1]
    for i in (firstindex(ps) + 2):lastindex(ps)
      !iscollinear(p₁, p₂, ps[i]) && return false
    end
  end
  return true
end

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

isconvex(::Frustum) = true

isconvex(::Torus) = false

isconvex(::Triangle) = true

isconvex(::Tetrahedron) = true

isconvex(p::Polygon) = _isconvex(p, Val(embeddim(p)))

_isconvex(p::Polygon, ::Val{2}) = Set(eachvertex(convexhull(p))) == Set(eachvertex(p))

_isconvex(p::Polygon, ::Val{3}) = isconvex(proj2D(p))

isconvex(h::Hexahedron) = all(isconvex, boundary(h))

isconvex(p::Pyramid) = isconvex(base(p))

isconvex(m::Multi) = isapproxequal(measure(convexhull(m)), measure(m))

# --------------
# OPTIMIZATIONS
# --------------

function isconvex(q::Quadrangle)
  A, B, C, D = vertices(q)
  d1 = Segment(A, C)
  d2 = Segment(B, D)
  intersects(d1, d2)
end
