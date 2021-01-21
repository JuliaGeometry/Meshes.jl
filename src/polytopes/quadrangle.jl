# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Quadrangle(p1, p2, p3, p4)

A quadrangle with points `p1`, `p2`, `p3`, `p4`.
"""
struct Quadrangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polygon{Dim,T}
  vertices::V
end

function measure(q::Quadrangle)
  vs = q.vertices
  Δ₁ = Triangle(view(vs, [1,2,3]))
  Δ₂ = Triangle(view(vs, [3,4,1]))
  measure(Δ₁) + measure(Δ₂)
end
