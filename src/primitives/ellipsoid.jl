# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(radii, center=(0, 0, 0), rotation=I)

A 3D ellipsoid with given `radii`, `center` and `rotation`.
"""
struct Ellipsoid{T,R} <: Primitive{3,T}
  radii::NTuple{3,T}
  center::Point{3,T}
  rotation::R
end

Ellipsoid(radii::NTuple{3,T}, center=(T(0), T(0), T(0)), rotation::R=I) where {T,R} =
  Ellipsoid{T,R}(radii, center, rotation)

paramdim(::Type{<:Ellipsoid}) = 2

radii(e::Ellipsoid) = e.radii

center(e::Ellipsoid) = e.center

rotation(e::Ellipsoid) = e.rotation

function (e::Ellipsoid{T})(θ, φ) where {T}
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

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ellipsoid{T}}) where {T} =
  Ellipsoid((rand(rng, T), rand(rng, T), rand(rng, T)), rand(rng, Point{3,T}), rand(rng, QuatRotation))
