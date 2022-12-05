# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(radius, bottom, top)

A solid circular cylinder embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    Cylinder(radius, segment)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` between origin of `bottom` and `top` planes.

    Cylinder(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
  radius::T
  bot::Plane{T}
  top::Plane{T}
end

function Cylinder(radius, segment::Segment{3,T}) where {T}
  a, b = extrema(segment)
  v    = b - a
  bot  = Plane(a, v)
  top  = Plane(b, v)
  Cylinder(T(radius), bot, top)
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

bottom(c::Cylinder) = c.bot

top(c::Cylinder) = c.top

center(c::Cylinder) = center(boundary(c))

axis(c::Cylinder) = axis(boundary(c))

isright(c::Cylinder) = isright(boundary(c))

boundary(c::Cylinder) = CylinderSurface(c.radius, c.bot, c.top)

measure(c::Cylinder) = norm(origin(c.bot) - origin(c.top)) * c.radius^2 * pi

volume(c::Cylinder) = measure(c)