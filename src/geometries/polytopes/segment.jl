# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment with end points `p1`, `p2`.
The segment can be called as `s(t)` with `t` between
`0` and `1` to interpolate linearly between its endpoints.

See also [`Rope`](@ref), [`Ring`](@ref), [`Line`](@ref).
"""
@polytope Segment 1 2

nvertices(::Type{<:Segment}) = 2

Base.minimum(s::Segment) = s.vertices[1]

Base.maximum(s::Segment) = s.vertices[2]

Base.extrema(s::Segment) = s.vertices[1], s.vertices[2]

function center(s::Segment)
  a, b = extrema(s)
  withcrs(s, (to(a) + to(b)) / 2)
end

==(s₁::Segment, s₂::Segment) = s₁.vertices == s₂.vertices

Base.isapprox(s₁::Segment, s₂::Segment; atol=atol(lentype(s₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(s₁.vertices, s₂.vertices))

function (s::Segment)(t)
  if t < 0 || t > 1
    throw(DomainError(t, "s(t) is not defined for t outside [0, 1]."))
  end
  a, b = s.vertices
  a + t * (b - a)
end

Base.reverse(s::Segment) = Segment(reverse(extrema(s)))
