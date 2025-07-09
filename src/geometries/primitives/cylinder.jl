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

See also [`CylinderSurface`](@ref).
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

bottom(c::Cylinder) = c.bot

top(c::Cylinder) = c.top

radius(c::Cylinder) = c.radius

function (c::Cylinder)(ρ, φ, z)
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((ρ, φ, z), "c(ρ, φ, z) is not defined for ρ, φ, z outside [0, 1]³."))
  end
  ℒ = lentype(c)
  T = numtype(ℒ)
  t = top(c)
  b = bottom(c)
  r = radius(c)
  a = axis(c)
  d = a(T(1)) - a(T(0))

  # rotation to align z axis with cylinder axis
  R = urotbetween(Vec(zero(ℒ), zero(ℒ), oneunit(ℒ)), d)

  # offset to translate cylinder to final position
  o = to(b(T(0), T(0)))

  # project a parametric segment between the top and bottom planes
  ρ′ = T(ρ) * r
  φ′ = T(φ) * 2 * T(π) * u"rad"
  p₁ = Point(convert(crs(c), Cylindrical(ρ′, φ′, zero(ℒ))))
  p₂ = Point(convert(crs(c), Cylindrical(ρ′, φ′, norm(d))))
  l = Line(p₁, p₂) |> Affine(R, o)
  s = Segment(l ∩ b, l ∩ t)
  s(T(z))
end

# ----------------------------------------------
# forward methods to boundary (CylinderSurface)
# ----------------------------------------------

axis(c::Cylinder) = axis(boundary(c))

isright(c::Cylinder) = isright(boundary(c))

hasintersectingplanes(c::Cylinder) = hasintersectingplanes(boundary(c))

==(c₁::Cylinder, c₂::Cylinder) = boundary(c₁) == boundary(c₂)

Base.isapprox(c₁::Cylinder, c₂::Cylinder; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(boundary(c₁), boundary(c₂); atol, kwargs...)
