# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    centroid(geometry)

The centroid of the `geometry`.
"""
centroid(g::Geometry) = center(g) # some geometries have a natural center

centroid(p::Point) = p

centroid(p::Polygon) = centroid(first(rings(p)))

centroid(p::Polytope) = withcrs(p, sum(to, vertices(p)) / nvertices(p))

centroid(b::Box) = withcrs(b, sum(to, extrema(b)) / 2)

centroid(p::Plane) = p(0, 0)

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
  x = vector.(1:n)
  w = volume.(1:n)
  all(iszero, w) && (w = ones(eltype(w), n))
  withcrs(d, sum(w .* x) / sum(w))
end

"""
    centroid(domain, ind)

The centroid of the `ind`-th element of the `domain`.
"""
centroid(d::Domain, ind::Int) = centroid(d[ind])

centroid(d::SubDomain, ind::Int) = centroid(parent(d), parentindices(d)[ind])

function centroid(g::QuasiCartesianGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  vertex(g, ijk) + Vec(spacing(g) ./ 2)
end

function centroid(g::RectilinearGrid{<:𝔼,<:CartesianOrProjected}, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p1 = vertex(g, ijk)
  p2 = vertex(g, ijk .+ 1)
  withcrs(g, (to(p1) + to(p2)) / 2)
end

centroid(m::TransformedMesh, ind::Int) = transform(m)(centroid(parent(m), ind))
