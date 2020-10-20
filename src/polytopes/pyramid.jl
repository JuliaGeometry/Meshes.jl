# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
struct Pyramid{Dim,T} <: Polytope{Dim,T,5}
    vertices::SVector{5,Point{Dim,T}}
end

