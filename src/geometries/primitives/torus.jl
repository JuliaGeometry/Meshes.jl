# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Torus(center, direction, major, minor)

A torus centered at `center` with axis of revolution directed by 
`direction` and with radii `major` and `minor`. 

"""
struct Torus{C<:CRS,Mₚ<:Manifold,V<:Vec{3},ℒ<:Len} <: Primitive{𝔼{3},C}
  center::Point{Mₚ,C}
  direction::V
  major::ℒ
  minor::ℒ
  Torus(center::Point{Mₚ,C}, direction::V, major::ℒ, minor::ℒ) where {C<:CRS,Mₚ<:Manifold,V<:Vec{3},ℒ<:Len} =
    new{C,Mₚ,V,float(ℒ)}(center, direction, major, minor)
end

Torus(center::Point, direction::Vec, major::Len, minor::Len) = Torus(center, direction, promote(major, minor)...)

Torus(center::Point, direction::Vec, major, minor) =
  Torus(center, direction, addunit(major, u"m"), addunit(minor, u"m"))

Torus(center::Tuple, direction::Tuple, major, minor) = Torus(Point(center), Vec(direction), major, minor)

"""
    Torus(p1, p2, p3, minor)

The torus whose centerline passes through points `p1`, `p2` and `p3` and with
minor radius `minor`.
"""
function Torus(p1::Point, p2::Point, p3::Point, minor::Len)
  c = Circle(p1, p2, p3)
  p = Plane(p1, p2, p3)
  Torus(center(c), normal(p), radius(c), minor)
end

Torus(p1::Point, p2::Point, p3::Point, minor) = Torus(p1, p2, p3, addunit(minor, u"m"))

Torus(p1::Tuple, p2::Tuple, p3::Tuple, minor) = Torus(Point(p1), Point(p2), Point(p3), minor)

paramdim(::Type{<:Torus}) = 2

center(t::Torus) = t.center

direction(t::Torus) = t.direction

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.direction)

==(t₁::Torus, t₂::Torus) =
  t₁.center == t₂.center && t₁.direction == t₂.direction && t₁.major == t₂.major && t₁.minor == t₂.minor

Base.isapprox(t₁::Torus, t₂::Torus; atol=atol(lentype(t₁)), kwargs...) =
  isapprox(t₁.center, t₂.center; atol, kwargs...) &&
  isapprox(t₁.direction, t₂.direction; atol, kwargs...) &&
  isapprox(t₁.major, t₂.major; atol, kwargs...) &&
  isapprox(t₁.minor, t₂.minor; atol, kwargs...)

function (t::Torus)(θ, φ)
  ℒ = lentype(t)
  T = numtype(ℒ)
  if (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((θ, φ), "t(θ, φ) is not defined for θ, φ outside [0, 1]²."))
  end
  c, n⃗ = t.center, t.direction
  R, r = t.major, t.minor

  Q = urotbetween(Vec(zero(ℒ), zero(ℒ), oneunit(ℒ)), n⃗)

  sθ, cθ = sincospi(2 * T(-θ))
  sφ, cφ = sincospi(2 * T(φ))
  x = (R + r * cθ) * cφ
  y = (R + r * cθ) * sφ
  z = r * sθ

  c + Q * Vec(x, y, z)
end
