# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    centroid(geometry)

The centroid of the `geometry`.
"""
centroid(g::Geometry) = withcrs(g, integral(to, g) / measure(g))

centroid(p::Point) = p

centroid(p::Plane) = center(p)

centroid(b::Box) = center(b)

centroid(b::Ball) = center(b)

centroid(s::Sphere) = center(s)

centroid(d::Disk) = center(d)

centroid(c::Circle) = center(c)

centroid(c::Cylinder) = centroid(boundary(c))

function centroid(c::CylinderSurface)
  a = centroid(bottom(c))
  b = centroid(top(c))
  withcrs(c, (to(a) + to(b)) / 2)
end

function centroid(p::ParaboloidSurface)
  c = apex(p)
  r = radius(p)
  f = focallength(p)
  z = r^2 / 4f
  x = zero(z)
  y = zero(z)
  c + Vec(x, y, z / 2)
end

centroid(t::Torus) = center(t)

centroid(s::Segment) = s(1 // 2)

centroid(t::Triangle) = t(1 // 3, 1 // 3)

centroid(q::Quadrangle) = q(1 // 2, 1 // 2)

centroid(t::Tetrahedron) = t(1 // 4, 1 // 4, 1 // 4)

centroid(h::Hexahedron) = h(1 // 2, 1 // 2, 1 // 2)

centroid(m::Multi) = centroid(GeometrySet(parent(m)))

centroid(g::TransformedGeometry) = transform(g)(centroid(parent(g)))

"""
    centroid(domain)

The centroid of the `domain`.
"""
function centroid(d::Domain)
  vector(i) = to(centroid(d, i))
  volume(i) = measure(element(d, i))
  n = nelements(d)
  x = map(vector, 1:n)
  w = map(volume, 1:n)
  all(iszero, w) && (w = ones(eltype(w), n))
  withcrs(d, sum(w .* x) / sum(w))
end

"""
    centroid(domain, ind)

The centroid of the `ind`-th element of the `domain`.
"""
centroid(d::Domain, ind::Int) = centroid(d[ind])

centroid(d::SubDomain, ind::Int) = centroid(parent(d), parentindices(d)[ind])

function centroid(g::OrthoRegularGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  vertex(g, ijk) + Vec(spacing(g) ./ 2)
end

function centroid(g::OrthoRectilinearGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = vertex(g, ijk)
  p2 = vertex(g, ijk .+ 1)
  withcrs(g, (to(p1) + to(p2)) / 2)
end

centroid(m::TransformedMesh, ind::Int) = transform(m)(centroid(parent(m), ind))
