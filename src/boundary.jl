# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    boundary(object)

Return the boundary of the `object`.
"""
function boundary end

boundary(::Point) = nothing

boundary(r::Ray) = r(0)

boundary(::Line) = nothing

function boundary(b::BezierCurve)
  p = controls(b)
  p₁, p₂ = first(p), last(p)
  p₁ ≈ p₂ ? nothing : Multi([p₁, p₂])
end

boundary(::Plane) = nothing

boundary(b::Box{1}) = Multi([minimum(b), maximum(b)])

function boundary(b::Box{2})
  A = coordinates(minimum(b))
  B = coordinates(maximum(b))
  v = Point.([(A[1], A[2]), (B[1], A[2]), (B[1], B[2]), (A[1], B[2])])
  Ring(v)
end

function boundary(b::Box{3})
  A = coordinates(minimum(b))
  B = coordinates(maximum(b))
  v =
    Point.([
      (A[1], A[2], A[3]),
      (B[1], A[2], A[3]),
      (B[1], B[2], A[3]),
      (A[1], B[2], A[3]),
      (A[1], A[2], B[3]),
      (B[1], A[2], B[3]),
      (B[1], B[2], B[3]),
      (A[1], B[2], B[3])
    ])
  c = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(v, connect.(c))
end

boundary(b::Ball) = Sphere(center(b), radius(b))

boundary(::Sphere) = nothing

boundary(::Ellipsoid) = nothing

boundary(d::Disk) = Circle(plane(d), radius(d))

boundary(::Circle) = nothing

boundary(c::Cylinder) = CylinderSurface(bottom(c), top(c), radius(c))

boundary(::CylinderSurface) = nothing

boundary(c::Cone) = ConeSurface(base(c), apex(c))

boundary(::ConeSurface) = nothing

boundary(f::Frustum) = FrustumSurface(bottom(f), top(f))

boundary(::FrustumSurface) = nothing

boundary(::ParaboloidSurface) = nothing

boundary(::Torus) = nothing

boundary(s::Segment) = Multi(pointify(s))

function boundary(r::Rope)
  v = vertices(r)
  Multi([first(v), last(v)])
end

boundary(::Ring) = nothing

boundary(p::Polygon) = hasholes(p) ? Multi(rings(p)) : first(rings(p))

function boundary(t::Tetrahedron)
  indices = [(3, 2, 1), (4, 1, 2), (4, 3, 1), (4, 2, 3)]
  SimpleMesh(pointify(t), connect.(indices))
end

function boundary(h::Hexahedron)
  indices = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(pointify(h), connect.(indices))
end

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(pointify(p), connect.(indices))
end

function boundary(m::Multi)
  bounds = [boundary(geom) for geom in parent(m)]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end
