# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ray(p, v)

A ray originating at point `p`, pointed in direction `v`.
It can be called as `r(t)` with `t > 0` to cast it at
`p + t * v`.
"""
struct Ray{Dim,T} <: Primitive{Dim,T}
  p::Point{Dim,T}
  v::Vec{Dim,T}
end

Ray(p::Tuple, v::Tuple) = Ray(Point(p), Vec(v))

paramdim(::Type{<:Ray}) = 1

isconvex(::Type{<:Ray}) = true

boundary(r::Ray) = r.p

function (r::Ray)(t)
  if t < 0
    throw(DomainError(t, "r(t) is not defined for t < 0."))
  end
  r.p + t * r.v
end

"""
    origin(ray)

The starting point of the ray.
"""
origin(r::Ray) = r.p

"""
    direction(ray)

The direction of the ray.
"""
direction(r::Ray) = r.v

function Base.in(p::Point{Dim, T}, r::Ray{Dim, T}) where {Dim,T}
	p ∈ Line(r(zero(T)), r(one(T))) && sum((p - r.p) .* r.v) >= 0 # additional check
end

==(r1::Ray{Dim, T}, r2::Ray{Dim, T}) where {Dim,T} = (r1.p == r2.p) && (r1(one(T)) ∈ r2) && (r2(one(T)) ∈ r1)
