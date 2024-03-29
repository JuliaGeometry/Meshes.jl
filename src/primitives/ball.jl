# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ball(center, radius)

A ball with `center` and `radius`.

See also [`Sphere`](@ref).
"""
struct Ball{Dim,T} <: Primitive{Dim,T}
  center::Point{Dim,T}
  radius::T
end

Ball(center::Point{Dim,T}, radius) where {Dim,T} = Ball(center, T(radius))

Ball(center::Tuple, radius) = Ball(Point(center), radius)

Ball(center::Point{Dim,T}) where {Dim,T} = Ball(center, T(1))

Ball(center::Tuple) = Ball(Point(center))

paramdim(::Type{<:Ball{Dim}}) where {Dim} = Dim

center(b::Ball) = b.center

radius(b::Ball) = b.radius

function (b::Ball{2,T})(ρ, φ) where {T}
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "b(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  c = b.center
  r = b.radius
  l = T(ρ) * r
  sφ, cφ = sincospi(2 * T(φ))
  x = l * cφ
  y = l * sφ
  c + Vec(x, y)
end

function (b::Ball{3,T})(ρ, θ, φ) where {T}
  if (ρ < 0 || ρ > 1) || (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, θ, φ), "b(ρ, θ, φ) is not defined for ρ, θ, φ outside [0, 1]³."))
  end
  c = b.center
  r = b.radius
  l = T(ρ) * r
  sθ, cθ = sincospi(T(θ))
  sφ, cφ = sincospi(2 * T(φ))
  x = l * sθ * cφ
  y = l * sθ * sφ
  z = l * cθ
  c + Vec(x, y, z)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ball{Dim,T}}) where {Dim,T} =
  Ball(rand(rng, Point{Dim,T}), rand(rng, T))
