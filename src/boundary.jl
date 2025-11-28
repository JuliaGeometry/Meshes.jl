# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    boundary(g)

The relative boundary of the geometry `g` is the subset of points
`p ∈ g` such that for all `ϵ > 0` the intersection of the `ϵ`-ball
centered at `p` with the affine hull of `g` contains points in `g`
*and* in the complement of `g`.

See <https://en.wikipedia.org/wiki/Relative_interior>.
"""
function boundary end

"""
    embedboundary(g)

The embedded boundary of the geometry `g` is the subset of points
`p ∈ g` such that for all `ϵ > 0` the `ϵ`-ball centered at `p`
contains points in `g` *and* in the complement of `g`.

See <https://en.wikipedia.org/wiki/Boundary_(topology)>.
"""
function embedboundary end

boundary(::Point) = nothing

embedboundary(p::Point) = p

boundary(r::Ray) = r(0)

embedboundary(r::Ray) = r

boundary(::Line) = nothing

embedboundary(l::Line) = l

function boundary(b::BezierCurve)
  p = controls(b)
  p₁, p₂ = first(p), last(p)
  p₁ ≈ p₂ ? nothing : Multi([p₁, p₂])
end

embedboundary(b::BezierCurve) = b

function boundary(c::ParametrizedCurve)
  p₁, p₂ = extrema(c)
  p₁ ≈ p₂ ? nothing : Multi([p₁, p₂])
end

embedboundary(c::ParametrizedCurve) = c

boundary(::Plane) = nothing

embedboundary(p::Plane) = p

boundary(b::Box{𝔼{1}}) = embedboundary(b)

boundary(b::Box{𝔼{2}}) = embedboundary(b)

boundary(b::Box{𝔼{3}}) = embedboundary(b)

boundary(b::Box{🌐}) = embedboundary(b)

embedboundary(b::Box{𝔼{1}}) = Multi([minimum(b), maximum(b)])

function embedboundary(b::Box{𝔼{2}})
  A = convert(Cartesian, coords(minimum(b)))
  B = convert(Cartesian, coords(maximum(b)))
  v = [withcrs(b, (A.x, A.y)), withcrs(b, (B.x, A.y)), withcrs(b, (B.x, B.y)), withcrs(b, (A.x, B.y))]
  Ring(v)
end

function embedboundary(b::Box{𝔼{3}})
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

function embedboundary(b::Box{🌐})
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

boundary(b::Ball) = embedboundary(b)

embedboundary(b::Ball) = Sphere(center(b), radius(b))

boundary(::Sphere) = nothing

embedboundary(s::Sphere) = s

boundary(::Ellipsoid) = nothing

embedboundary(e::Ellipsoid) = e

boundary(d::Disk) = Circle(plane(d), radius(d))

embedboundary(d::Disk) = d

boundary(::Circle) = nothing

embedboundary(c::Circle) = c

boundary(c::Cylinder) = embedboundary(c)

embedboundary(c::Cylinder) = CylinderSurface(bottom(c), top(c), radius(c))

boundary(::CylinderSurface) = nothing

embedboundary(c::CylinderSurface) = c

boundary(c::Cone) = embedboundary(c)

embedboundary(c::Cone) = ConeSurface(base(c), apex(c))

boundary(::ConeSurface) = nothing

embedboundary(c::ConeSurface) = c

boundary(f::Frustum) = embedboundary(f)

embedboundary(f::Frustum) = FrustumSurface(bottom(f), top(f))

boundary(::FrustumSurface) = nothing

embedboundary(f::FrustumSurface) = f

boundary(::ParaboloidSurface) = nothing

embedboundary(p::ParaboloidSurface) = p

boundary(::Torus) = nothing

embedboundary(t::Torus) = t

boundary(s::Segment) = Multi([minimum(s), maximum(s)])

embedboundary(s::Segment) = paramdim(s) < embeddim(s) ? s : boundary(s)

function boundary(r::Rope)
  v = vertices(r)
  Multi([first(v), last(v)])
end

embedboundary(r::Rope) = paramdim(r) < embeddim(r) ? r : boundary(r)

boundary(::Ring) = nothing

embedboundary(r::Ring) = paramdim(r) < embeddim(r) ? r : boundary(r)

boundary(p::Polygon) = hasholes(p) ? Multi(rings(p)) : first(rings(p))

embedboundary(p::Polygon) = paramdim(p) < embeddim(p) ? p : boundary(p)

boundary(t::Tetrahedron) = embedboundary(t)

function embedboundary(t::Tetrahedron)
  indices = [(3, 2, 1), (4, 1, 2), (4, 3, 1), (4, 2, 3)]
  SimpleMesh(collect(eachvertex(t)), connect.(indices))
end

boundary(h::Hexahedron) = embedboundary(h)

function embedboundary(h::Hexahedron)
  indices = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(collect(eachvertex(h)), connect.(indices))
end

boundary(p::Pyramid) = embedboundary(p)

function embedboundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(collect(eachvertex(p)), connect.(indices))
end

boundary(w::Wedge) = embedboundary(w)

function embedboundary(w::Wedge)
  indices = [(1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)]
  SimpleMesh(collect(eachvertex(w)), connect.(indices))
end

function boundary(m::Multi)
  bounds = [boundary(geom) for geom in parent(m)]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

function embedboundary(m::Multi)
  bounds = [embedboundary(geom) for geom in parent(m)]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

function boundary(g::TransformedGeometry)
  b = boundary(parent(g))
  isnothing(b) ? nothing : transform(g)(b)
end

function embedboundary(g::TransformedGeometry)
  b = embedboundary(parent(g))
  isnothing(b) ? nothing : transform(g)(b)
end
