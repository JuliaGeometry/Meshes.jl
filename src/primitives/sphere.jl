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

center(s::Sphere) = s.center
radius(s::Sphere) = s.radius

function Base.in(p::Point, s::Sphere)
  x = coordinates(p)
  c = coordinates(s.center)
  r = s.radius
  sum(abs2, x - c) ≈ r^2
end
