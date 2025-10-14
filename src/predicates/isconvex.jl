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

isconvex(h::Hexahedron) = _isconvex(h)

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

function _isconvex(h::Hexahedron)
  b = boundary(h)
  all(isconvex, b) && _isgloballyconvex(b)
end

function _isgloballyconvex(b)
  t = convert(HalfEdgeTopology, topology(b))
  # map segments to their adjacent faces
  𝒞₁₂ = Coboundary{1,2}(t)
  r = 0
  for edgeᵢ in 1:nfacets(t)
    adjfaces = 𝒞₁₂(edgeᵢ)
    length(adjfaces) == 2 || return false
    # normal faces
    n₁, n₂ = map(adjfaces) do faceᵢ
      face = element(b, faceᵢ)
      verts = vertices(face)
      n = length(verts)
      n ≥ 3 ? normal(Plane(verts[1], verts[2], verts[3])) : throw(ArgumentError("Face must have at least 3 vertices"))
    end

    # calculate dihedral angle
    # for convex polyhedra, all but 2 should be positive (I believe)
    cosθ = (n₁ ⋅ n₂) / (norm(n₁) * norm(n₂))
    r += cosθ < -1e-10
    r > 2 && return false  # early exit if too many seam edges
  end
  true
end
