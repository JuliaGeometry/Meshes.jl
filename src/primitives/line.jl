# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(a, b)

A line passing through points `a` and `b`.
"""
struct Line{Dim,T} <: Primitive{Dim,T}
  a::Point{Dim,T}
  b::Point{Dim,T}
end

paramdim(::Type{<:Line}) = 1
