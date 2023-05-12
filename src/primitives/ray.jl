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

measure(::Ray{Dim,T}) where {Dim,T} = typemax(T)

Base.length(r::Ray) = measure(r)

boundary(r::Ray) = r.p

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

Base.in(p::Point, r::Ray) = p ∈ Line(r.p, r.p + r.v) && (p - r.p) ⋅ r.v ≥ 0

==(r₁::Ray, r₂::Ray) = (r₁.p ≈ r₂.p) && (r₁.p + r₁.v) ∈ r₂ && (r₂.p + r₂.v) ∈ r₁

function (r::Ray)(t)
  if t < 0
    throw(DomainError(t, "r(t) is not defined for t < 0."))
  end
  r.p + t * r.v
end
