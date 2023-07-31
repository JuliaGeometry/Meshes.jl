# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(a, b)

A line passing through points `a` and `b`.

See also [`Segment`](@ref).
"""
struct Line{Dim,T} <: Primitive{Dim,T}
  a::Point{Dim,T}
  b::Point{Dim,T}
end

Line(a::Tuple, b::Tuple) = Line(Point(a), Point(b))

paramdim(::Type{<:Line}) = 1

measure(::Line{Dim,T}) where {Dim,T} = typemax(T)

Base.length(l::Line) = measure(l)

perimeter(::Line{Dim,T}) where {Dim,T} = zero(T)

boundary(::Line) = nothing

==(l1::Line, l2::Line) = l1.a ∈ l2 && l1.b ∈ l2 && l2.a ∈ l1 && l2.b ∈ l1

(l::Line)(t) = l.a + t * (l.b - l.a)

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Line{Dim,T}}) where {Dim,T} =
  Line(rand(rng, Point{Dim,T}, 2)...)
