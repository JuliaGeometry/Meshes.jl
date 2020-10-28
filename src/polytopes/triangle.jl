# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Triangle(p1, p2, p3)

A triangle with points `p1`, `p2`, `p3`.
"""
struct Triangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{Dim,T}
  vertices::V
end

function volume(triangle::Triangle)
  A, B, C = triangle.vertices
  abs((B - A) × (C - A)) / 2
end

function Base.in(P::Point, triangle::Triangle)
  A, B, C = triangle.vertices
  bc = C - B
  ca = A - C
  ab = B - A
  ap = P - A
  bp = P - B
  cp = P - C

  abp = bc[1] * bp[2] - bc[2] * bp[1]
  cap = ab[1] * ap[2] - ab[2] * ap[1]
  bcp = ca[1] * cp[2] - ca[2] * cp[1]

  (abp ≥ 0) && (bcp ≥ 0) && (cap ≥ 0)
end
