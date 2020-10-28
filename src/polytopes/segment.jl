# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{Dim,T}
  vertices::V
end
