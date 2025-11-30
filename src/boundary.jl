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

embedboundary(b::Box{𝔼{1}}) = boundary(b)

embedboundary(b::Box{𝔼{2}}) = boundary(b)

embedboundary(b::Box{𝔼{3}}) = boundary(b)

embedboundary(b::Box{🌐}) = boundary(b)

boundary(b::Ball) = Sphere(center(b), radius(b))

embedboundary(b::Ball) = boundary(b)

boundary(::Sphere) = nothing

embedboundary(s::Sphere) = s

boundary(::Ellipsoid) = nothing

embedboundary(e::Ellipsoid) = e

boundary(d::Disk) = Circle(plane(d), radius(d))

embedboundary(d::Disk) = d

boundary(::Circle) = nothing

embedboundary(c::Circle) = c

boundary(c::Cylinder) = CylinderSurface(bottom(c), top(c), radius(c))

embedboundary(c::Cylinder) = boundary(c)

boundary(::CylinderSurface) = nothing

embedboundary(c::CylinderSurface) = c

boundary(c::Cone) = ConeSurface(base(c), apex(c))

embedboundary(c::Cone) = boundary(c)

boundary(::ConeSurface) = nothing

embedboundary(c::ConeSurface) = c

boundary(f::Frustum) = FrustumSurface(bottom(f), top(f))

embedboundary(f::Frustum) = boundary(f)

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

function boundary(t::Tetrahedron)
  indices = [(3, 2, 1), (4, 1, 2), (4, 3, 1), (4, 2, 3)]
  SimpleMesh(collect(eachvertex(t)), connect.(indices))
end

embedboundary(t::Tetrahedron) = boundary(t)

function boundary(h::Hexahedron)
  indices = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(collect(eachvertex(h)), connect.(indices))
end

embedboundary(h::Hexahedron) = boundary(h)

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(collect(eachvertex(p)), connect.(indices))
end

embedboundary(p::Pyramid) = boundary(p)

function boundary(w::Wedge)
  indices = [(1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)]
  SimpleMesh(collect(eachvertex(w)), connect.(indices))
end

embedboundary(w::Wedge) = boundary(w)

function boundary(m::Multi)
  bs = filter(!isnothing, [boundary(g) for g in parent(m)])
  isempty(bs) ? nothing : reduce(merge, bs)
end

function embedboundary(m::Multi)
  bs = [embedboundary(g) for g in parent(m)]
  reduce(merge, bs)
end

function boundary(g::TransformedGeometry)
  b = boundary(parent(g))
  t = transform(g)
  if isnothing(b)
    nothing
  elseif b isa Geometry
    TransformedGeometry(b, t)
  elseif b isa Mesh
    TransformedMesh(b, t)
  end
end

function embedboundary(g::TransformedGeometry)
  b = embedboundary(parent(g))
  t = transform(g)
  if b isa Geometry
    TransformedGeometry(b, t)
  elseif b isa Mesh
    TransformedMesh(b, t)
  end
end

"""
    boundarypoints(geometry)

Return vector of [`Point`](@ref)s that approximate
the [`embedboundary`](@ref) of the `geometry`.
"""
boundarypoints(g::Geometry) = _boundarypoints(embedboundary(g))

# discretize boundary and extract vertices
_boundarypoints(g::Geometry) = vertices(discretize(g))

# --------------
# OPTIMIZATIONS
# --------------

_boundarypoints(p::Point) = [p]

_boundarypoints(m::MultiPoint) = parent(m)

_boundarypoints(p::Polytope) = _boundarypoints(manifold(p), p)

_boundarypoints(::Type{<:𝔼}, p::Polytope) = vertices(p)

_boundarypoints(::Type{<:🌐}, p::Polytope) = vertices(discretize(p))

_boundarypoints(m::Mesh) = vertices(m)
