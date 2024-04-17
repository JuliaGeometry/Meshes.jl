# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Plane(p, u, v)

A plane embedded in R³ passing through point `p`,
defined by non-parallel vectors `u` and `v`.

    Plane(p, n)

Alternatively specify point `p` and a given normal
vector `n` to the plane.
"""
struct Plane{P<:Point,T} <: Primitive
  p::P
  u::Vec{3,T}
  v::Vec{3,T}
end

function Plane{T}(p::P, n::Vec{3,T}) where {P<:Point,T}
  u, v = householderbasis(n)
  Plane{P,T}(p, u, v)
end

Plane(p::Point, n::Vec{3,T}) where {T} = Plane{T}(p, n)

Plane(p::Tuple, u::Tuple, v::Tuple) = Plane(Point(p), Vec(u), Vec(v))

Plane(p::Tuple, n::Tuple) = Plane(Point(p), Vec(n))

# TODO
# function Plane(p1::Point{3,T}, p2::Point{3,T}, p3::Point{3,T}) where {T}
#   t = Triangle(p1, p2, p3)
#   if isapprox(area(t), zero(T), atol=atol(T))
#     throw(ArgumentError("The three points are colinear."))
#   end
#   Plane(p1, normal(t))
# end

paramdim(::Type{<:Plane}) = 2

normal(p::Plane) = normalize(p.u × p.v)

==(p₁::Plane, p₂::Plane) =
  p₁(0, 0) ∈ p₂ && p₁(1, 0) ∈ p₂ && p₁(0, 1) ∈ p₂ && p₂(0, 0) ∈ p₁ && p₂(1, 0) ∈ p₁ && p₂(0, 1) ∈ p₁

Base.isapprox(p₁::Plane{P,T}, p₂::Plane{P,T}) where {P,T} =
  isapprox((p₁(0, 0) - p₂(0, 0)) ⋅ normal(p₂), zero(T), atol=atol(T)) &&
  isapprox((p₂(0, 0) - p₁(0, 0)) ⋅ normal(p₁), zero(T), atol=atol(T)) &&
  isapprox(_area(normal(p₁), normal(p₂)), zero(T), atol=atol(T))

_area(v₁::Vec, v₂::Vec) = norm(v₁ × v₂)

(p::Plane)(u, v) = p.p + u * p.u + v * p.v

# TODO
# Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Plane{T}}) where {T} =
#   Plane(rand(rng, Point{3,T}), rand(rng, Vec{3,T}))
