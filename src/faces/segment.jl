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
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T,1}
  vertices::V
end
