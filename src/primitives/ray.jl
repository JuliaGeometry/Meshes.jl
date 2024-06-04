# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ray(p, v)

A ray originating at point `p`, pointed in direction `v`.
It can be called as `r(t)` with `t > 0` to cast it at
`p + t * v`.
"""
struct Ray{Dim,C<:CRS,ℒ<:Len} <: Primitive{Dim,C}
  p::Point{Dim,C}
  v::Vec{Dim,ℒ}
end

Ray(p::Tuple, v::Tuple) = Ray(Point(p), Vec(v))

paramdim(::Type{<:Ray}) = 1

==(r₁::Ray, r₂::Ray) = (r₁.p ≈ r₂.p) && (r₁.p + r₁.v) ∈ r₂ && (r₂.p + r₂.v) ∈ r₁

function (r::Ray)(t)
  if t < 0
    throw(DomainError(t, "r(t) is not defined for t < 0."))
  end
  r.p + t * r.v
end

Random.rand(rng::Random.AbstractRNG, ::Type{Ray{Dim}}) where {Dim} =
  Ray(rand(rng, Point{Dim}), rand(rng, Vec{Dim,Met{Float64}}))
