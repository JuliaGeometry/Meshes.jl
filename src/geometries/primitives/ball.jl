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

Ball(center::Point, radius) = Ball(center, aslen(radius))

Ball(center::Tuple, radius) = Ball(Point(center), radius)

Ball(center::Point) = Ball(center, oneunit(lentype(center)))

Ball(center::Tuple) = Ball(Point(center))

paramdim(::Type{<:Ball{𝔼{Dim}}}) where {Dim} = Dim

paramdim(::Type{<:Ball{🌐}}) = 2

center(b::Ball) = b.center

radius(b::Ball) = b.radius

==(b₁::Ball, b₂::Ball) = center(b₁) == center(b₂) && radius(b₁) == radius(b₂)

Base.isapprox(b₁::Ball, b₂::Ball; atol=atol(lentype(b₁)), kwargs...) =
  isapprox(center(b₁), center(b₂); atol, kwargs...) && isapprox(radius(b₁), radius(b₂); atol, kwargs...)

function (b::Ball{𝔼{2}})(ρ, φ)
  C = basecrs(b)
  T = numtype(lentype(b))
  ρ′ = ρ * radius(b)
  φ′ = φ * 2 * T(π) * u"rad"
  p = Point(convert(C, Polar(ρ′, φ′)))
  p + to(center(b))
end

function (b::Ball{𝔼{3}})(ρ, θ, φ)
  C = basecrs(b)
  T = numtype(lentype(b))
  ρ′ = ρ * radius(b)
  θ′ = θ * T(π) * u"rad"
  φ′ = φ * 2 * T(π) * u"rad"
  p = Point(convert(C, Spherical(ρ′, θ′, φ′)))
  p + to(center(b))
end
