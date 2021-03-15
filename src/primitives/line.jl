# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(p1, p2)

A line passing through points `p1` and `p2`.

See also [`Segment`](@ref).
"""
struct Line{Dim,T} <: Primitive{Dim,T}
  a::Point{Dim,T}
  b::Point{Dim,T}
end

paramdim(::Type{<:Line}) = 1
