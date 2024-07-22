# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ray(p, v)

A ray originating at point `p`, pointed in direction `v`.
It can be called as `r(t)` with `t > 0` to cast it at
`p + t * v`.
"""
struct Ray{M<:𝔼,C<:CRS,V<:Vec} <: Primitive{M,C}
  p::Point{M,C}
  v::V
end

Ray(p::Tuple, v::Tuple) = Ray(Point(p), Vec(v))

paramdim(::Type{<:Ray}) = 1

==(r₁::Ray, r₂::Ray) = r₁.p == r₂.p && (r₁.p + r₁.v) ∈ r₂ && (r₂.p + r₂.v) ∈ r₁

Base.isapprox(r₁::Ray, r₂::Ray; atol=atol(lentype(r₁)), kwargs...) =
  isapprox(r₁.p, r₂.p; atol, kwargs...) && isapprox(r₁.v, r₂.v; atol, kwargs...)

function (r::Ray)(t)
  if t < 0
    throw(DomainError(t, "r(t) is not defined for t < 0."))
  end
  r.p + t * r.v
end
