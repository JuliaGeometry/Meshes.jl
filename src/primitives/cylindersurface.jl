# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderSurface(radius, axis, bottom, top)

A circular cylinder surface embedded in R³ with given `radius` around `axis`,
delimited by `bottom` and `top` planes.

    CylinderSurface(radius, segment)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` of axis.

    CylinderSurface(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct CylinderSurface{T} <: Primitive{3,T}
  radius::T
  axis::Line{3,T}
  bot::Plane{T}
  top::Plane{T}
end

function CylinderSurface(radius::T, segment::Segment{3,T}) where {T}
  a, b = extrema(segment)
  v    = b - a
  axis = Line(a, b)
  bot  = Plane(a, v)
  top  = Plane(b, v)
  CylinderSurface(radius, axis, bot, top)
end

function CylinderSurface(radius::T) where {T}
  _0 = (zero(T), zero(T), zero(T))
  _1 = (zero(T), zero(T), one(T))
  segment = Segment(_0, _1)
  CylinderSurface(radius, segment)
end

paramdim(::Type{<:CylinderSurface}) = 2

isconvex(::Type{<:CylinderSurface}) = true

radius(c::CylinderSurface) = c.radius

axis(c::CylinderSurface) = c.axis

boundary(::CylinderSurface) = nothing