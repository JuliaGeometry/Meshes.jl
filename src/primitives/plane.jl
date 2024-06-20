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
struct Plane{C<:CRS,ℒ<:Len} <: Primitive{3,C}
  p::Point{3,C}
  u::Vec{3,ℒ}
  v::Vec{3,ℒ}
end

function Plane(p::Point{3}, n::Vec{3})
  u, v = householderbasis(n)
  Plane(p, u, v)
end

Plane(p::Tuple, u::Tuple, v::Tuple) = Plane(Point(p), Vec(u), Vec(v))

Plane(p::Tuple, n::Tuple) = Plane(Point(p), Vec(n))

function Plane(p1::Point{3}, p2::Point{3}, p3::Point{3})
  t = Triangle(p1, p2, p3)
  if isapproxzero(area(t))
    throw(ArgumentError("The three points are colinear."))
  end
  Plane(p1, normal(t))
end

paramdim(::Type{<:Plane}) = 2

normal(p::Plane) = unormalize(ucross(p.u, p.v))

==(p₁::Plane, p₂::Plane) =
  p₁(0, 0) ∈ p₂ && p₁(1, 0) ∈ p₂ && p₁(0, 1) ∈ p₂ && p₂(0, 0) ∈ p₁ && p₂(1, 0) ∈ p₁ && p₂(0, 1) ∈ p₁

Base.isapprox(p₁::Plane, p₂::Plane; atol=atol(lentype(p₁)), kwargs...) =
  isapproxzero(udot(p₁(0, 0) - p₂(0, 0), normal(p₂)); atol, kwargs...) &&
  isapproxzero(udot(p₂(0, 0) - p₁(0, 0), normal(p₁)); atol, kwargs...) &&
  isapproxzero(norm(ucross(normal(p₁), normal(p₂))); atol, kwargs...)

(p::Plane)(u, v) = p.p + u * p.u + v * p.v

Random.rand(rng::Random.AbstractRNG, ::Type{Plane}) = Plane(rand(rng, Point{3}), rand(rng, Vec{3,Met{Float64}}))
