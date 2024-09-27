# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(a, b)

A line passing through points `a` and `b`.

See also [`Segment`](@ref).
"""
struct Line{M<:Manifold,C<:CRS} <: Primitive{M,C}
  a::Point{M,C}
  b::Point{M,C}
end

Line(a::Tuple, b::Tuple) = Line(Point(a), Point(b))

paramdim(::Type{<:Line}) = 1

==(l₁::Line, l₂::Line) = l₁.a ∈ l₂ && l₁.b ∈ l₂ && l₂.a ∈ l₁ && l₂.b ∈ l₁

Base.isapprox(l₁::Line, l₂::Line; atol=atol(lentype(l₁)), kwargs...) =
  isapprox(l₁.a, l₂.a; atol, kwargs...) && isapprox(l₁.b, l₂.b; atol, kwargs...)

(l::Line)(t) = coordsum((l.a, l.b), weights=((1 - t), t))
