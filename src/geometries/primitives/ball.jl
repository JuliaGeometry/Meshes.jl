# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ball(center, radius)

A ball with `center` and `radius`.

See also [`Sphere`](@ref).
"""
struct Ball{M<:Manifold,C<:CRS,ℒ<:Len} <: Primitive{M,C}
  center::Point{M,C}
  radius::ℒ
  Ball(center::Point{M,C}, radius::ℒ) where {M<:Manifold,C<:CRS,ℒ<:Len} = new{M,C,float(ℒ)}(center, radius)
end

Ball(center::Point, radius) = Ball(center, addunit(radius, u"m"))

Ball(center::Tuple, radius) = Ball(Point(center), radius)

Ball(center::Point) = Ball(center, oneunit(lentype(center)))

Ball(center::Tuple) = Ball(Point(center))

paramdim(::Type{<:Ball{𝔼{Dim}}}) where {Dim} = Dim

paramdim(::Type{<:Ball{🌐}}) = 2

center(b::Ball) = b.center

radius(b::Ball) = b.radius

==(b₁::Ball, b₂::Ball) = b₁.center == b₂.center && b₁.radius == b₂.radius

Base.isapprox(b₁::Ball, b₂::Ball; atol=atol(lentype(b₁)), kwargs...) =
  isapprox(b₁.center, b₂.center; atol, kwargs...) && isapprox(b₁.radius, b₂.radius; atol, kwargs...)

(b::Ball)(args...) = _ball(Val(embeddim(b)), b, args...)

function _ball(::Val{2}, b, ρ, φ)
  T = numtype(lentype(b))
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

function _ball(::Val{3}, b, ρ, θ, φ)
  T = numtype(lentype(b))
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
