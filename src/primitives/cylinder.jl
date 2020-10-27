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

radius(c::Cylinder) = c.radius
height(c::Cylinder) = norm(c.finish - c.start)
volume(c::Cylinder) = π * radius(c)^2 * height(c)
