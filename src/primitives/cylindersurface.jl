# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderSurface(radius, bottom, top)

A circular cylinder surface embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    CylinderSurface(radius, segment)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` between origin of `bottom` and `top` planes.

    CylinderSurface(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct CylinderSurface{T} <: Primitive{3,T}
  radius::T
  bot::Plane{T}
  top::Plane{T}
end

function CylinderSurface(radius::T, segment::Segment{3,T}) where {T}
  a, b = extrema(segment)
  v    = b - a
  bot  = Plane(a, v)
  top  = Plane(b, v)
  CylinderSurface(radius, bot, top)
end

function CylinderSurface(radius::T) where {T}
  _0 = (T(0), T(0), T(0))
  _1 = (T(0), T(0), T(1))
  segment = Segment(_0, _1)
  CylinderSurface(radius, segment)
end

paramdim(::Type{<:CylinderSurface}) = 2

isconvex(::Type{<:CylinderSurface}) = true

radius(c::CylinderSurface) = c.radius

axis(c::CylinderSurface) = Line(origin(c.bot), origin(c.top))

planes(c::CylinderSurface) = (c.bot, c.top)

function isright(c::CylinderSurface{T}) where {T}
  # cylinder is right if axis
  # is aligned with plane normals
  a = axis(c)
  d = a(T(1)) - a(T(0))
  v = normal(c.bot)
  w = normal(c.top)
  isparallelv = isapprox(norm(d × v), zero(T), atol=atol(T))
  isparallelw = isapprox(norm(d × w), zero(T), atol=atol(T))
  isparallelv && isparallelw
end

boundary(::CylinderSurface) = nothing