# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Torus(center, normal, major, minor)

A torus centered at `center` with axis of revolution directed by 
`normal` and with radii `major` and `minor`. 

"""
struct Torus{T} <: Primitive{3,T}
  center::Point{3,T}
  normal::Vec{3,T}
  major::T
  minor::T
end

Torus(center::Point{3,T}, normal::Vec{3,T}, major, minor) where {T} = Torus(center, normal, T(major), T(minor))

Torus(center::Tuple, normal::Tuple, major, minor) = Torus(Point(center), Vec(normal), major, minor)

"""
    Torus(p1, p2, p3, minor)

The torus whose equator passes through points `p1`, `p2` and `p3` and with
minor radius `minor`.
"""
function Torus(p1::Point{3,T}, p2::Point{3,T}, p3::Point{3,T}, minor) where {T}
  c = Circle(p1, p2, p3)
  p = Plane(p1, p2, p3)
  Torus(center(c), normal(p), radius(c), T(minor))
end

Torus(p1::Tuple, p2::Tuple, p3::Tuple, minor) = Torus(Point(p1), Point(p2), Point(p3), minor)

paramdim(::Type{<:Torus}) = 2

isconvex(::Type{<:Torus}) = false

isperiodic(::Type{<:Torus}) = (true, true)

center(t::Torus) = t.center

normal(t::Torus) = t.normal

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.normal)

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  R, r = t.major, t.minor
  4π^2 * R * r
end

area(t::Torus) = measure(t)

function Base.in(p::Point, t::Torus)
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor
  M = uvrotation(n⃗, Vec(0, 0, 1))
  x, y, z = M * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end
