# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(radii, center=(0, 0, 0), rotation=I)

A 3D ellipsoid with given `radii`, `center` and `rotation`.
"""
struct Ellipsoid{ℒ<:Len,C<:CRS,R} <: Primitive{3,C}
  radii::NTuple{3,ℒ}
  center::Point{3,C}
  rotation::R
  Ellipsoid{ℒ,C,R}(radii, center, rotation) where {ℒ<:Len,C<:CRS,R} = new(radii, center, rotation)
end

Ellipsoid(radii::NTuple{3,ℒ}, center::Point{3,C}, rotation::R) where {ℒ<:Len,C<:CRS,R} =
  Ellipsoid{float(ℒ),C,R}(radii, center, rotation)

Ellipsoid(radii::NTuple{3}, center::Point{3}, rotation) = Ellipsoid(addunit.(radii, u"m"), center, rotation)

Ellipsoid(radii::NTuple{3}, center::NTuple{3}, rotation) = Ellipsoid(radii, Point(center), rotation)

Ellipsoid(radii::NTuple{3,T}, center=(zero(T), zero(T), zero(T)), rotation=I) where {T} =
  Ellipsoid(radii, center, rotation)

paramdim(::Type{<:Ellipsoid}) = 2

radii(e::Ellipsoid) = e.radii

center(e::Ellipsoid) = e.center

rotation(e::Ellipsoid) = e.rotation

Base.isapprox(e₁::Ellipsoid, e₂::Ellipsoid) =
  all(e₁.radii .≈ e₂.radii) && e₁.center ≈ e₂.center && e₁.rotation ≈ e₂.rotation

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

Random.rand(rng::Random.AbstractRNG, ::Type{Ellipsoid}) = Ellipsoid(
  (rand(rng, Met{Float64}), rand(rng, Met{Float64}), rand(rng, Met{Float64})),
  rand(rng, Point{3}),
  rand(rng, QuatRotation)
)
