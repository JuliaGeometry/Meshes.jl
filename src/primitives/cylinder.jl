# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(radius, axis, bottom, top)

A solid circular cylinder embedded in R³ with given `radius` around `axis`,
delimited by `bottom` and `top` planes.

    Cylinder(radius, segment)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` of axis.

    Cylinder(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
  radius::T
  axis::Line{3,T}
  bot::Plane{T}
  top::Plane{T}
end

function Cylinder(radius::T, segment::Segment{3,T}) where {T}
  a, b = extrema(segment)
  v    = b - a
  axis = Line(a, b)
  bot  = Plane(a, v)
  top  = Plane(b, v)
  Cylinder(radius, axis, bot, top)
end

function Cylinder(radius::T) where {T}
  _0 = (T(0), T(0), T(0))
  _1 = (T(0), T(0), T(1))
  segment = Segment(_0, _1)
  Cylinder(radius, segment)
end

paramdim(::Type{<:Cylinder}) = 3

isconvex(::Type{<:Cylinder}) = true

radius(c::Cylinder) = c.radius

axis(c::Cylinder) = c.axis

planes(c::Cylinder) = (c.bot, c.top)

isright(c::Cylinder) = isright(boundary(c))

boundary(c::Cylinder) =
  CylinderSurface(c.radius, c.axis, c.bot, c.top)