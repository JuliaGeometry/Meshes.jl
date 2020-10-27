# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Quadrangle(p1, p2, p3, p4)

A quadrangle with points `p1`, `p2`, `p3`, `p4`.
"""
struct Quadrangle{Dim,T} <: Polytope{Dim,T,4}
  vertices::NTuple{4,Point{Dim,T}}
end
