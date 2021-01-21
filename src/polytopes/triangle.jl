# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Triangle(p1, p2, p3)

A triangle with points `p1`, `p2`, `p3`.
"""
struct Triangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polygon{Dim,T}
  vertices::V
end

"""
    signarea(triangle)

Compute the signed area of `triangle`.
"""
function signarea(t::Triangle{2})
  a, b, c = coordinates.(t.vertices)
  ((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])) / 2
end

measure(t::Triangle{2}) = abs(signarea(t))
