# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(radii, center=(0, 0, 0), rotation=I)

A 3D ellipsoid with given `radii`, `center` and `rotation`.
"""
struct Ellipsoid{ℒ<:Len,P<:Point{3},R} <: Primitive{3}
  radii::NTuple{3,ℒ}
  center::P
  rotation::R
  Ellipsoid{ℒ,P,R}(radii, center, rotation) where {ℒ<:Len,P<:Point{3},R} = new(radii, center, rotation)
end

Ellipsoid(radii::NTuple{3,ℒ}, center::P, rotation::R) where {ℒ<:Len,P<:Point{3},R} =
  Ellipsoid{float(ℒ),P,R}(radii, center, rotation)

Ellipsoid(radii::NTuple{3}, center::P, rotation::R) where {P<:Point{3},R} =
  Ellipsoid(addunit.(radii, u"m"), center, rotation)

Ellipsoid(radii::NTuple{3,T}, center=(zero(T), zero(T), zero(T)), rotation=I) where {T} =
  Ellipsoid(radii, Point(center), rotation)

paramdim(::Type{<:Ellipsoid}) = 2

lentype(::Type{<:Ellipsoid{ℒ,P}}) where {ℒ,P} = lentype(P)

radii(e::Ellipsoid) = e.radii

center(e::Ellipsoid) = e.center

rotation(e::Ellipsoid) = e.rotation

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
