# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ray(p, v)

A ray originating at point `p`, pointed in direction `v`.
It can be called as `r(t)` with `t > 0` to cast it at
`p + t * v`.
"""
struct Ray{C<:CRS,Mₚ<:Manifold,Dim,V<:Vec{Dim}} <: Primitive{𝔼{Dim},C}
  p::Point{Mₚ,C}
  v::V
end

Ray(p::Tuple, v::Tuple) = Ray(Point(p), Vec(v))

paramdim(::Type{<:Ray}) = 1

==(r₁::Ray, r₂::Ray) = r₁.p == r₂.p && (r₁.p + r₁.v) ∈ r₂ && (r₂.p + r₂.v) ∈ r₁

Base.isapprox(r₁::Ray, r₂::Ray; atol=atol(lentype(r₁)), kwargs...) =
  isapprox(r₁.p, r₂.p; atol, kwargs...) && isapprox(r₁.v, r₂.v; atol, kwargs...)

(r::Ray)(t) = r.p + t * r.v
