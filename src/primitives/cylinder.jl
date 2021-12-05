# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(radius[, axis, bottom, top])

A solid circular cylinder with given `radius` around `axis`
that is delimited with `bottom` and `top` planes.

## Parameters

* `radius` - Radius of circular base
* `axis`   - Axis of orientation (default to `Line((0,0,0), (0,0,1))`)
* `bottom` - Bottom plane (default to `Plane((0,0,0), (0,0,1))`)
* `top`    - Top plane (default to `Plane((0,0,1), (0,0,1))`)

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
  radius::T
  axis::Line{3,T}
  bot::Plane{T}
  top::Plane{T}
end

function Cylinder(radius::T) where {T}
  _0   = (zero(T), zero(T), zero(T))
  _1   = (zero(T), zero(T), one(T))
  axis = Line( _0, _1)
  bot  = Plane(_0, _1)
  top  = Plane(_1, _1)
  Cylinder(radius, axis, bot, top)
end

paramdim(::Type{<:Cylinder}) = 3

isconvex(::Type{<:Cylinder}) = true

radius(c::Cylinder) = c.radius