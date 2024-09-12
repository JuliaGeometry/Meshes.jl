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
struct Cylinder{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{𝔼{3},C}
  bot::P
  top::P
  radius::ℒ
  Cylinder(bot::P, top::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new{C,P,float(ℒ)}(bot, top, radius)
end

Cylinder(bot::P, top::P, radius) where {P<:Plane} = Cylinder(bot, top, addunit(radius, u"m"))

function Cylinder(start::Point, finish::Point, radius)
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  Cylinder(bot, top, radius)
end

Cylinder(start::Tuple, finish::Tuple, radius) = Cylinder(Point(start), Point(finish), radius)

Cylinder(start::Point, finish::Point) = Cylinder(start, finish, oneunit(lentype(start)))

Cylinder(start::Tuple, finish::Tuple) = Cylinder(Point(start), Point(finish))

function Cylinder(radius)
  z = zero(radius)
  o = oneunit(radius)
  Cylinder(Point(z, z, z), Point(z, z, o), radius)
end

paramdim(::Type{<:Cylinder}) = 3

radius(c::Cylinder) = c.radius

bottom(c::Cylinder) = c.bot

top(c::Cylinder) = c.top

axis(c::Cylinder) = axis(boundary(c))

isright(c::Cylinder) = isright(boundary(c))

hasintersectingplanes(c::Cylinder) = hasintersectingplanes(boundary(c))

==(c₁::Cylinder, c₂::Cylinder) = boundary(c₁) == boundary(c₂)

Base.isapprox(c₁::Cylinder, c₂::Cylinder; kwargs...) = isapprox(boundary(c₁), boundary(c₂); kwargs...)

function (c::Cylinder)(ρ, φ, z)
  ℒ = lentype(c)
  T = numtype(ℒ)
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((ρ, φ, z), "c(ρ, φ, z) is not defined for ρ, φ, z outside [0, 1]³."))
  end
  t = top(c)
  b = bottom(c)
  r = radius(c)
  a = axis(c)
  d = a(T(1)) - a(T(0))
  h = norm(d)
  o = b(0, 0)

  # rotation to align z axis with cylinder axis
  Q = urotbetween(Vec(zero(ℒ), zero(ℒ), oneunit(ℒ)), d)

  # project a parametric segment between the top and bottom planes
  lsφ, lcφ = T(ρ) * r .* sincospi(2 * T(φ))
  p₁ = o + Q * Vec(lcφ, lsφ, zero(ℒ))
  p₂ = o + Q * Vec(lcφ, lsφ, h)
  l = Line(p₁, p₂)
  s = Segment(l ∩ b, l ∩ t)
  s(T(z))
end
