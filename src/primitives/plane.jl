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

function Plane(p1::Point{3,T1}, p2::Point{3,T2}, p3::Point{3,T3}) where {T1,T2,T3}
  T = promote_type(T1, T2, T3)
  n⃗ = Vec{3,T}((p2 - p1) × (p3 - p1))
  if isapprox(area(Triangle(p1, p2, p3)), zero(T), atol = atol(T))
    throw(ArgumentError("The three points are colinear."))
  end
  Plane(convert(Point{3,T}, p1), n⃗)
end

paramdim(::Type{<:Plane}) = 2

isconvex(::Type{<:Plane}) = true

measure(::Plane{T}) where {T} = typemax(T)

area(p::Plane) = measure(p)

boundary(::Plane) = nothing

perimeter(::Plane{T}) where {T} = zero(T)

(p::Plane)(s, t) = p.p + s*p.u + t*p.v

"""
    origin(plane)

Return the origin of the `plane`.
"""
origin(p::Plane) = p.p

"""
    normal(plane)

Normal vector to the `plane`.
"""
normal(p::Plane) = normalize(p.u × p.v)

==(p₁::Plane{T}, p₂::Plane{T}) where {T} =
  isapprox((p₁.v - p₁.u) ⋅ normal(p₂), zero(T), atol=atol(T)) &&
  isapprox((p₂.v - p₂.u) ⋅ normal(p₁), zero(T), atol=atol(T))

function Base.in(pt::Point, pl::Plane)
  d = normal(pl) ⋅ (pt - origin(pl))
  T = typeof(d)
  isapprox(d, zero(T), atol = atol(T))
end
  