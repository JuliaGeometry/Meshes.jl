# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderSurface(radius, bottom, top)

A circular cylinder surface embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    CylinderSurface(radius, segment)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` between `bottom` and `top` planes.

    CylinderSurface(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct CylinderSurface{T} <: Primitive{3,T}
  radius::T
  bot::Plane{T}
  top::Plane{T}
end

function CylinderSurface(radius, segment::Segment{3,T}) where {T}
  a, b = extrema(segment)
  v = b - a
  bot = Plane(a, v)
  top = Plane(b, v)
  return CylinderSurface(T(radius), bot, top)
end

function CylinderSurface(radius::T) where {T}
  _0 = (T(0), T(0), T(0))
  _1 = (T(0), T(0), T(1))
  segment = Segment(_0, _1)
  return CylinderSurface(radius, segment)
end

paramdim(::Type{<:CylinderSurface}) = 2

isconvex(::Type{<:CylinderSurface}) = true

radius(c::CylinderSurface) = c.radius

bottom(c::CylinderSurface) = c.bot

top(c::CylinderSurface) = c.top

function center(c::CylinderSurface)
  a = coordinates(c.bot(0, 0))
  b = coordinates(c.top(0, 0))
  return Point((a .+ b) ./ 2)
end

axis(c::CylinderSurface) = Line(c.bot(0, 0), c.top(0, 0))

function isright(c::CylinderSurface{T}) where {T}
  # cylinder is right if axis
  # is aligned with plane normals
  a = axis(c)
  d = a(T(1)) - a(T(0))
  v = normal(c.bot)
  w = normal(c.top)
  isparallelv = isapprox(norm(d × v), zero(T); atol=atol(T))
  isparallelw = isapprox(norm(d × w), zero(T); atol=atol(T))
  return isparallelv && isparallelw
end

boundary(::CylinderSurface) = nothing

function measure(c::CylinderSurface{T}) where {T}
  return (norm(c.bot(0, 0) - c.top(0, 0)) + c.radius) * 2 * c.radius * T(π)
end

area(c::CylinderSurface) = measure(c)
