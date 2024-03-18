# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(bottom, top, radius)

A solid circular cylinder embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    Cylinder(start, finish, radius)

Alternatively, construct a right circular cylinder with given `radius`
along the segment with `start` and `finish` end points.

    Cylinder(start, finish)

Or construct a right circular cylinder with unit radius along the segment
with `start` and `finish` end points.

    Cylinder(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See <https://en.wikipedia.org/wiki/Cylinder>. 
"""
struct Cylinder{T} <: Primitive{3,T}
  bot::Plane{T}
  top::Plane{T}
  radius::T
end

function Cylinder(start::Point{3,T}, finish::Point{3,T}, radius) where {T}
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  Cylinder(bot, top, T(radius))
end

Cylinder(start::Tuple, finish::Tuple, radius) = Cylinder(Point(start), Point(finish), radius)

Cylinder(start::Point{3,T}, finish::Point{3,T}) where {T} = Cylinder(start, finish, T(1))

Cylinder(start::Tuple, finish::Tuple) = Cylinder(Point(start), Point(finish))

Cylinder(radius::T) where {T} = Cylinder(Point(T(0), T(0), T(0)), Point(T(0), T(0), T(1)), radius)

paramdim(::Type{<:Cylinder}) = 3

radius(c::Cylinder) = c.radius

bottom(c::Cylinder) = c.bot

top(c::Cylinder) = c.top

center(c::Cylinder) = center(boundary(c))

axis(c::Cylinder) = axis(boundary(c))

isright(c::Cylinder) = isright(boundary(c))

hasintersectingplanes(c::Cylinder) = hasintersectingplanes(boundary(c))

Base.isapprox(c₁::Cylinder, c₂::Cylinder) = boundary(c₁) ≈ boundary(c₂)

function (c::Cylinder{T})(ρ, φ, z) where {T}
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((ρ, φ, z), "c(ρ, φ, z) is not defined for ρ, φ, z outside [0, 1]³."))
  end

  t = top(c)
  b = bottom(c)
  r = radius(c)
  a = axis(c)
  d = a(T(1)) - a(T(0))
  h = norm(d)

  # calculate translation/rotation to map between cylinder-space and global coords
  o = b(0, 0)
  Q = rotation_between(Vec{3,T}(0, 0, 1), d)

  # project a parametric Segment between the top and bottom planes
  lsφ, lcφ = T(ρ) * r .* sincospi(2 * T(φ))
  p₁ = o + Q * Vec(lcφ, lsφ, T(0))
  p₂ = o + Q * Vec(lcφ, lsφ, h)
  l = Line(p₁, p₂)
  s = Segment(l ∩ b, l ∩ t)
  s(T(z))
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Cylinder{T}}) where {T} =
  Cylinder(rand(rng, Plane{T}), rand(rng, Plane{T}), rand(rng, T))
