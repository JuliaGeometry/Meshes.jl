# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment from point `p1` to `p2`.

See also [`Geodesic`](@ref), [`Line`](@ref).
"""
@polytope Segment 1 2

nvertices(::Type{<:Segment}) = 2

Base.minimum(s::Segment) = s.vertices[1]

Base.maximum(s::Segment) = s.vertices[2]

Base.extrema(s::Segment) = s.vertices[1], s.vertices[2]

center(s::Segment) = coordmean(extrema(s))

==(s₁::Segment, s₂::Segment) = s₁.vertices == s₂.vertices

Base.isapprox(s₁::Segment, s₂::Segment; atol=atol(lentype(s₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(s₁.vertices, s₂.vertices))

function (s::Segment)(t)
  a, b = s.vertices
  coordsum((a, b), weights=((1 - t), t))
end

Base.reverse(s::Segment) = Segment(reverse(extrema(s)))
