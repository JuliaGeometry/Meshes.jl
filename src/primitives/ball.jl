# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ball(center, radius)

A ball with `center` and `radius`.

See also [`Sphere`](@ref).
"""
struct Ball{Dim,P<:Point{Dim},L<:Len} <: Primitive{Dim}
  center::P
  radius::L
  Ball(center::P, radius::L) where {Dim,P<:Point{Dim},L<:Len} = new{Dim,P,float(L)}(center, radius)
end

Ball(center::Point, radius) = Ball(center, addunit(radius, u"m"))

Ball(center::Tuple, radius) = Ball(Point(center), radius)

Ball(center::Point) = Ball(center, 1.0u"m")

Ball(center::Tuple) = Ball(Point(center))

paramdim(::Type{<:Ball{Dim}}) where {Dim} = Dim

center(b::Ball) = b.center

radius(b::Ball) = b.radius

function (b::Ball{2,P,L})(ρ, φ) where {P,L}
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "b(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  c = b.center
  r = b.radius
  l = L(ρ) * r
  sφ, cφ = sincospi(2 * L(φ))
  x = l * cφ
  y = l * sφ
  c + Vec(x, y)
end

function (b::Ball{3,P,L})(ρ, θ, φ) where {P,L}
  if (ρ < 0 || ρ > 1) || (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, θ, φ), "b(ρ, θ, φ) is not defined for ρ, θ, φ outside [0, 1]³."))
  end
  c = b.center
  r = b.radius
  l = L(ρ) * r
  sθ, cθ = sincospi(L(θ))
  sφ, cφ = sincospi(2 * L(φ))
  x = l * sθ * cφ
  y = l * sθ * sφ
  z = l * cθ
  c + Vec(x, y, z)
end

# TODO
# Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ball{Dim,T}}) where {Dim,T} =
#   Ball(rand(rng, Point{Dim,T}), rand(rng, T))
