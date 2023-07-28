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
struct Plane{T} <: Primitive{3,T}
  p::Point{3,T}
  u::Vec{3,T}
  v::Vec{3,T}
end

function Plane{T}(p::Point{3,T}, n::Vec{3,T}) where {T}
  u, v = householderbasis(n)
  Plane{T}(p, u, v)
end

Plane(p::Point{3,T}, n::Vec{3,T}) where {T} = Plane{T}(p, n)

Plane(p::Tuple, u::Tuple, v::Tuple) = Plane(Point(p), Vec(u), Vec(v))

Plane(p::Tuple, n::Tuple) = Plane(Point(p), Vec(n))

function Plane(p1::Point{3,T}, p2::Point{3,T}, p3::Point{3,T}) where {T}
  t = Triangle(p1, p2, p3)
  if isapprox(area(t), zero(T), atol=atol(T))
    throw(ArgumentError("The three points are colinear."))
  end
  Plane(p1, normal(t))
end

paramdim(::Type{<:Plane}) = 2

isconvex(::Type{<:Plane}) = true

isparametrized(::Type{<:Plane}) = true

normal(p::Plane) = normalize(p.u × p.v)

measure(::Plane{T}) where {T} = typemax(T)

area(p::Plane) = measure(p)

perimeter(::Plane{T}) where {T} = zero(T)

boundary(::Plane) = nothing

Base.isapprox(p₁::Plane{T}, p₂::Plane{T}) where {T} =
  isapprox((p₁.v - p₁.u) ⋅ normal(p₂), zero(T), atol=atol(T)) &&
  isapprox((p₂.v - p₂.u) ⋅ normal(p₁), zero(T), atol=atol(T))

Base.in(pt::Point{3,T}, pl::Plane{T}) where {T} = isapprox(normal(pl) ⋅ (pt - pl(0, 0)), zero(T), atol=atol(T))

(p::Plane)(u, v) = p.p + u * p.u + v * p.v

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Plane{T}}) where {T} =
  Plane(rand(rng, Point{3,T}), rand(rng, Vec{3,T}))
