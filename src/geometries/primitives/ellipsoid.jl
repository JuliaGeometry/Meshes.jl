# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(radii, center=(0, 0, 0), rotation=I)

A 3D ellipsoid with given `radii`, `center` and `rotation`.
"""
struct Ellipsoid{C<:CRS,Mₚ<:Manifold,R,ℒ<:Len} <: Primitive{𝔼{3},C}
  radii::NTuple{3,ℒ}
  center::Point{Mₚ,C}
  rotation::R
  Ellipsoid(radii::NTuple{3,ℒ}, center::Point{Mₚ,C}, rotation::R) where {C<:CRS,Mₚ<:Manifold,R,ℒ<:Len} =
    new{C,Mₚ,R,float(ℒ)}(radii, center, rotation)
end

Ellipsoid(radii::Tuple, center::Point, rotation) = Ellipsoid(addunit.(radii, u"m"), center, rotation)

Ellipsoid(radii::Tuple, center::Tuple, rotation) = Ellipsoid(radii, Point(center), rotation)

Ellipsoid(radii::Tuple, center=(_zero(radii), _zero(radii), _zero(radii)), rotation=I) =
  Ellipsoid(radii, center, rotation)

_zero(radii) = zero(first(radii))

paramdim(::Type{<:Ellipsoid}) = 2

radii(e::Ellipsoid) = e.radii

center(e::Ellipsoid) = e.center

rotation(e::Ellipsoid) = e.rotation

==(e₁::Ellipsoid, e₂::Ellipsoid) = radii(e₁) == radii(e₂) && center(e₁) == center(e₂) && rotation(e₁) == rotation(e₂)

function Base.isapprox(e₁::Ellipsoid, e₂::Ellipsoid; atol=atol(lentype(e₁)), kwargs...)
  u = Unitful.promote_unit(unit(lentype(e₁)), unit(lentype(e₂)))
  all(isapprox(r₁, r₂; atol, kwargs...) for (r₁, r₂) in zip(radii(e₁), radii(e₂))) &&
    isapprox(center(e₁), center(e₂); atol, kwargs...) &&
    isapprox(rotation(e₁), rotation(e₂); atol=ustrip(u, atol), kwargs...)
end

function (e::Ellipsoid)(θ, φ)
  if (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((θ, φ), "e(θ, φ) is not defined for θ, φ outside [0, 1]²."))
  end
  T = numtype(lentype(e))
  r = radii(e)
  c = center(e)
  R = rotation(e)
  sθ, cθ = sincospi(T(θ))
  sφ, cφ = sincospi(2 * T(φ))
  x = r[1] * sθ * cφ
  y = r[2] * sθ * sφ
  z = r[3] * cθ
  c + R * Vec(x, y, z)
end
