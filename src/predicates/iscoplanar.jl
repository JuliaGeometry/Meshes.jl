# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    iscoplanar(A, B, C, D)

Tells whether or not the points `A`, `B`, `C` and `D` are coplanar.
"""
iscoplanar(A::Point{3}, B::Point{3}, C::Point{3}, D::Point{3}) = isapproxzero(volume(Tetrahedron(A, B, C, D)))
