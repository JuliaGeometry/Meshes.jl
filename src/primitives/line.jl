# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(a, b)

A line passing through points `a` and `b`.

See also [`Segment`](@ref).
"""
struct Line{Dim,P<:Point} <: Primitive{Dim}
  a::P
  b::P
end

Line(a::Tuple, b::Tuple) = Line(Point(a), Point(b))

paramdim(::Type{<:Line}) = 1

coordtype(::Type{<:Line{Dim,P}}) where {Dim,P} = coordtype(P)

==(l₁::Line, l₂::Line) = l₁.a ∈ l₂ && l₁.b ∈ l₂ && l₂.a ∈ l₁ && l₂.b ∈ l₁

(l::Line)(t) = l.a + t * (l.b - l.a)

# TODO
# Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Line{Dim,T}}) where {Dim,T} =
#   Line(rand(rng, Point{Dim,T}, 2)...)
