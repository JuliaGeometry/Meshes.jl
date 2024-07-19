# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(radii, center=(0, 0, 0), rotation=I)

A 3D ellipsoid with given `radii`, `center` and `rotation`.
"""
struct Ellipsoid{M<:AbstractManifold,C<:CRS,R,ℒ<:Len} <: Primitive{M,C}
  radii::NTuple{3,ℒ}
  center::Point{M,C}
  rotation::R
  Ellipsoid(radii::NTuple{3,ℒ}, center::Point{M,C}, rotation::R) where {M<:AbstractManifold,C<:CRS,R,ℒ<:Len} =
    new{M,C,R,float(ℒ)}(radii, center, rotation)
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

==(e₁::Ellipsoid, e₂::Ellipsoid) = e₁.radii == e₂.radii && e₁.center == e₂.center && e₁.rotation == e₂.rotation

function Base.isapprox(e₁::Ellipsoid, e₂::Ellipsoid; atol=atol(lentype(e₁)), kwargs...)
  u = Unitful.promote_unit(unit(lentype(e₁)), unit(lentype(e₂)))
  all(isapprox(r₁, r₂; atol, kwargs...) for (r₁, r₂) in zip(e₁.radii, e₂.radii)) &&
    isapprox(e₁.center, e₂.center; atol, kwargs...) &&
    isapprox(e₁.rotation, e₂.rotation; atol=ustrip(u, atol), kwargs...)
end

function (e::Ellipsoid)(θ, φ)
  T = numtype(lentype(e))
  if (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((θ, φ), "e(θ, φ) is not defined for θ, φ outside [0, 1]²."))
  end
  r = e.radii
  c = e.center
  R = e.rotation
  sθ, cθ = sincospi(T(θ))
  sφ, cφ = sincospi(2 * T(φ))
  x = r[1] * sθ * cφ
  y = r[2] * sθ * sφ
  z = r[3] * cθ
  c + R * Vec(x, y, z)
end
