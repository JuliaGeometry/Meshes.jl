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
measure(t::Triangle{3}) = area(t)

function Base.in(p::Point{2}, t::Triangle{2})
  a, b, c = t.vertices
  abp = signarea(a, b, p)
  bcp = signarea(b, c, p)
  cap = signarea(c, a, p)
  areas = (abp, bcp, cap)
  all(areas .≥ 0) || all(areas .≤ 0)
end

function Base.in(p::Point{3}, t::Triangle{3})
  # https://people.cs.clemson.edu/~dhouse/courses/404/notes/barycentric.pdf
  a, b, c = t.vertices

  v_n = (b - a) × (c - b)

  A = norm(v_n)

  n = v_n / A

  u = ((c - b) × (p - b)) ⋅ n / A
  v = ((a - c) × (p - c)) ⋅ n / A
  
  (u ≥ 0) && (v ≥ 0) && ((u + v) ≤ 1)
end