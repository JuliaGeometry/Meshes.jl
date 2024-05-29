# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(a, b)

A line passing through points `a` and `b`.

See also [`Segment`](@ref).
"""
struct Line{Dim,P<:Point{Dim}} <: Primitive{Dim}
  a::P
  b::P
end

Line(a::Tuple, b::Tuple) = Line(Point(a), Point(b))

paramdim(::Type{<:Line}) = 1

lentype(::Type{<:Line{Dim,P}}) where {Dim,P} = lentype(P)

==(l₁::Line, l₂::Line) = l₁.a ∈ l₂ && l₁.b ∈ l₂ && l₂.a ∈ l₁ && l₂.b ∈ l₁

(l::Line)(t) = l.a + t * (l.b - l.a)

Random.rand(rng::Random.AbstractRNG, ::Type{Line{Dim}}) where {Dim} = Line(rand(rng, Point{Dim}), rand(rng, Point{Dim}))
