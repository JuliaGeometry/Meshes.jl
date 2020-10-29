# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ball(center, radius)

A Ball with `center` and `radius`.
"""
struct Ball{Dim,T} <: Primitive{Dim,T}
  center::Point{Dim,T}
  radius::T
end

Ball(center::Tuple, radius) = Ball(Point(center), radius)

center(b::Ball) = b.center
radius(b::Ball) = b.radius

function Base.in(p::Point, b::Ball)
  x = coordinates(p)
  c = coordinates(b.center)
  r = b.radius
  sum(abs2, x - c) ≤ r^2
end
