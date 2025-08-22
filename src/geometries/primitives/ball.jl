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

function (b::Ball{𝔼{2}})(ρ, φ)
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "b(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  T = numtype(lentype(b))
  ρ′ = T(ρ) * radius(b)
  φ′ = T(φ) * 2 * T(π) * u"rad"
  p = Point(convert(crs(b), Polar(ρ′, φ′)))
  p + to(center(b))
end

function (b::Ball{𝔼{3}})(ρ, θ, φ)
  if (ρ < 0 || ρ > 1) || (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, θ, φ), "b(ρ, θ, φ) is not defined for ρ, θ, φ outside [0, 1]³."))
  end
  T = numtype(lentype(b))
  ρ′ = T(ρ) * radius(b)
  θ′ = T(θ) * T(π) * u"rad"
  φ′ = T(φ) * 2 * T(π) * u"rad"
  p = Point(convert(crs(b), Spherical(ρ′, θ′, φ′)))
  p + to(center(b))
end
