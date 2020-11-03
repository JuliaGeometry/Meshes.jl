# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
#       v
#       ^
#       |
#       |
# 0-----+-----1 --> u

"""
    Segment(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T}
  vertices::V
end

paramdim(::Type{<:Segment}) = 1

facets(s::Segment) = (v for v in s.vertices)
