# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Hexahedron(p1, p2, ..., p8)

A hexahedron with points `p1`, `p2`, ..., `p8`.
"""
struct Hexahedron{Dim,T} <: Polytope{Dim,T,8}
  vertices::NTuple{8,Point{Dim,T}}
end
