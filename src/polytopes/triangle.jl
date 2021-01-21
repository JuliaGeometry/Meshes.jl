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

measure(t::Triangle{2}) = abs(signarea(t))

function Base.in(p::Point{2}, t::Triangle{2})
  a, b, c = t.vertices
  abp = signarea(a, b, p)
  bcp = signarea(b, c, p)
  cap = signarea(c, a, p)
  areas = (abp, bcp, cap)
  all(areas .≥ 0) || all(areas .≤ 0)
end

"""
    signarea(triangle)

Compute the signed area of `triangle`.
"""
signarea(t::Triangle{2}) = signarea(t.vertices[1], t.vertices[2], t.vertices[3])

"""
    signarea(p₁, p₂, p₃)

Compute signed area of triangle formed
by points `p₁`, `p₂` and `p₃`.
"""
function signarea(p₁::Point{2}, p₂::Point{2}, p₃::Point{2})
  a = coordinates(p₁)
  b = coordinates(p₂)
  c = coordinates(p₃)
  ((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])) / 2
end
