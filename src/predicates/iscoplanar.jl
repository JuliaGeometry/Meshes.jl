# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    iscoplanar(A, B, C, D)

Tells whether or not the points `A`, `B`, `C` and `D` are coplanar.
"""
iscoplanar(A::Point, B::Point, C::Point, D::Point) = isapproxzero(volume(Tetrahedron(A, B, C, D)))
