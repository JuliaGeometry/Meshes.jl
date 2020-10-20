# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Line(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Line{Dim,T} <: Polytope{Dim,T,2}
    vertices::SVector{2,Point{Dim,T}}
end
