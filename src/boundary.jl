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

boundary(b::Box{𝔼{1}}) = Multi([minimum(b), maximum(b)])

function boundary(b::Box{𝔼{2}})
  A = convert(Cartesian, coords(minimum(b)))
  B = convert(Cartesian, coords(maximum(b)))
  v = [withcrs(b, (A.x, A.y)), withcrs(b, (B.x, A.y)), withcrs(b, (B.x, B.y)), withcrs(b, (A.x, B.y))]
  Ring(v)
end

function boundary(b::Box{𝔼{3}})
  A = convert(Cartesian, coords(minimum(b)))
  B = convert(Cartesian, coords(maximum(b)))
  v = [
    withcrs(b, (A.x, A.y, A.z)),
    withcrs(b, (B.x, A.y, A.z)),
    withcrs(b, (B.x, B.y, A.z)),
    withcrs(b, (A.x, B.y, A.z)),
    withcrs(b, (A.x, A.y, B.z)),
    withcrs(b, (B.x, A.y, B.z)),
    withcrs(b, (B.x, B.y, B.z)),
    withcrs(b, (A.x, B.y, B.z))
  ]
  c = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(v, connect.(c))
end

function boundary(b::Box{🌐})
  A = convert(LatLon, coords(minimum(b)))
  B = convert(LatLon, coords(maximum(b)))
  v = [
    withcrs(b, (A.lat, A.lon), LatLon),
    withcrs(b, (A.lat, B.lon), LatLon),
    withcrs(b, (B.lat, B.lon), LatLon),
    withcrs(b, (B.lat, A.lon), LatLon)
  ]
  Ring(v)
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

function boundary(w::Wedge)
  indices = [(1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)]
  SimpleMesh(pointify(w), connect.(indices))
end

function boundary(m::Multi)
  bounds = [boundary(geom) for geom in parent(m)]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

boundary(g::TransformedGeometry) = transform(g)(boundary(parent(g)))
