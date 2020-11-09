# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Triangle(p1, p2, p3)

A triangle with points `p1`, `p2`, `p3`.
"""
struct Triangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{2,Dim,T}
  vertices::V
end