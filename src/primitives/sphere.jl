# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sphere(center, radius)

A sphere with `center` and `radius`.
"""
struct Sphere{Dim,T} <: Primitive{Dim,T}
  center::Point{Dim,T}
  radius::T
end

Sphere(center::Tuple, radius) = Sphere(Point(center), radius)

"""
    Sphere(p1, p2, p3)

A 2D sphere passing through points `p1`, `p2` and `p3`.
"""
function Sphere(p1::Point{2}, p2::Point{2}, p3::Point{2})
  x1, y1 = p2 - p1
  x2, y2 = p3 - p2
  c1 = centroid(Segment(p1, p2))
  c2 = centroid(Segment(p2, p3))
  l1 = Line(c1, c1 + Vec(y1, -x1))
  l2 = Line(c2, c2 + Vec(y2, -x2))
  center = l1 ∩ l2
  radius = norm(center - p2)
  Sphere(center, radius)
end

Sphere(p1::Tuple, p2::Tuple, p3::Tuple) = Sphere(Point(p1), Point(p2), Point(p3))

paramdim(::Type{<:Sphere{Dim}}) where {Dim} = Dim - 1

isconvex(::Type{<:Sphere}) = false

center(s::Sphere) = s.center
radius(s::Sphere) = s.radius

# https://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
function measure(s::Sphere{Dim}) where {Dim}
  r, n = s.radius, Dim
  2π^(n/2)*r^(n-1) / gamma(n/2)
end

function Base.in(p::Point, s::Sphere)
  x = coordinates(p)
  c = coordinates(s.center)
  r = s.radius
  sum(abs2, x - c) ≈ r^2
end
