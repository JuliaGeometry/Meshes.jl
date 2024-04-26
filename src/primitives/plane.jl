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
struct Plane{P<:Point{3},V<:Vec{3}} <: Primitive{3}
  p::P
  u::V
  v::V
end

function Plane(p::Point{3}, n::Vec{3})
  u, v = householderbasis(n)
  Plane(p, u, v)
end

Plane(p::Tuple, u::Tuple, v::Tuple) = Plane(Point(p), Vec(u), Vec(v))

Plane(p::Tuple, n::Tuple) = Plane(Point(p), Vec(n))

function Plane(p1::Point{3}, p2::Point{3}, p3::Point{3})
  t = Triangle(p1, p2, p3)
  a = area(t)
  if isapprox(a, zero(a), atol=atol(a))
    throw(ArgumentError("The three points are colinear."))
  end
  Plane(p1, normal(t))
end

paramdim(::Type{<:Plane}) = 2

lentype(::Type{<:Plane{P}}) where {P} = lentype(P)

normal(p::Plane) = Vec(normalize(p.u × p.v) * unit(lentype(p)))

==(p₁::Plane, p₂::Plane) =
  p₁(0, 0) ∈ p₂ && p₁(1, 0) ∈ p₂ && p₁(0, 1) ∈ p₂ && p₂(0, 0) ∈ p₁ && p₂(1, 0) ∈ p₁ && p₂(0, 1) ∈ p₁

function Base.isapprox(p₁::Plane, p₂::Plane)
  x₁ = (p₁(0, 0) - p₂(0, 0)) ⋅ normal(p₂)
  if isapprox(x₁, zero(x₁), atol=atol(x₁))
    x₂ = (p₂(0, 0) - p₁(0, 0)) ⋅ normal(p₁)
    if isapprox(x₂, zero(x₂), atol=atol(x₂))
      x₃ = _area(normal(p₁), normal(p₂))
      if isapprox(x₃, zero(x₃), atol=atol(x₃))
        return true
      end
    end
  end
  false
end

_area(v₁::Vec, v₂::Vec) = norm(v₁ × v₂)

(p::Plane)(u, v) = p.p + u * p.u + v * p.v

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Plane}) =
  Plane(rand(rng, Point{3}), rand(rng, Vec{3,Met{Float64}}))
