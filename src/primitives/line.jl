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

isconvex(::Type{<:Line}) = true

(l::Line)(t) = l.a + t * (l.b - l.a)

function Base.in(p::Point, l::Line)
  w = norm(l.b-l.a)
  d = evaluate(Euclidean(), p, l)
  # d ≈ 0.0 will be too precise, and d < atol{T} can't scale.
  d+w ≈ w
end

==(l1::Line, l2::Line) = l1.a ∈ l2 && l1.b ∈ l2 && l2.a ∈ l1 && l2.b ∈ l1
