# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(start, finish, radius)

A right circular cylinder with `start` and `finish` points,
and `radius` of revolution. See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
  start::Point{3,T}
  finish::Point{3,T}
  radius::T
end

Cylinder(start::Tuple, finish::Tuple, radius) =
  Cylinder(Point(start), Point(finish), radius)

paramdim(::Type{<:Cylinder}) = 3

isconvex(::Type{<:Cylinder}) = true

radius(c::Cylinder) = c.radius
height(c::Cylinder) = norm(c.finish - c.start)

measure(c::Cylinder) = π * radius(c)^2 * height(c)
