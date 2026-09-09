# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment from point `p1` to point `p2`.

See also [`Rope`](@ref), [`Ring`](@ref), [`Line`](@ref).
"""
struct Segment{M<:Manifold,C<:CRS} <: Chain{M,C}
  a::Point{M,C}
  b::Point{M,C}
end

Segment(a::Tuple, b::Tuple) = Segment(Point(a), Point(b))

Segment(ab::AbstractVector) = Segment(ab...)

Segment(ab::Tuple) = Segment(ab...)

nvertices(::Type{<:Segment}) = 2

vertices(s::Segment) = SVector(s.a, s.b)

Base.minimum(s::Segment) = s.a

Base.maximum(s::Segment) = s.b

Base.extrema(s::Segment) = s.a, s.b

center(s::Segment{<:𝔼}) = s(1 // 2)

==(s₁::Segment, s₂::Segment) = s₁.a == s₂.a && s₁.b == s₂.b

Base.isapprox(s₁::Segment, s₂::Segment; atol=atol(lentype(s₁)), kwargs...) =
  isapprox(s₁.a, s₂.a; atol=atol, kwargs...) && isapprox(s₁.b, s₂.b; atol=atol, kwargs...)

(s::Segment{<:𝔼})(t) = s.a + t * (s.b - s.a)

function (s::Segment{🌐})(t)
  d = GeodesicDistance()
  ϕ = geodesicbwd(s.a, s.b)
  geodesicfwd(s.a, ϕ, t * d(s.a, s.b))
end

Base.reverse(s::Segment) = Segment(s.b, s.a)
