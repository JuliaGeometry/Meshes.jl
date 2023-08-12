# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    iscoplanar(A, B, C, D)

Tells whether or not the points `A`, `B`, `C` and `D` are coplanar.
"""
function iscoplanar(A::Point{3,T}, B::Point{3,T}, C::Point{3,T}, D::Point{3,T}) where {T}
  vol = volume(Tetrahedron(A, B, C, D))
  isapprox(vol, zero(T), atol=atol(T))
end
