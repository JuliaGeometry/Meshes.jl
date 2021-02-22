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

function measure(tri::Triangle{3})
	v1, v2, v3 = vertices(tri)
	e1 = v2 - v1
	e2 = v3 - v1
	c = cross(e1, e2)
	norm(c) * 0.5
end

function normal(tri::Triangle{3})
    v1, v2, v3 = vertices(tri)
	e1 = v2 - v1
	e2 = v3 - v1
	Vec(normalize(cross(e1, e2)))
end

function Base.in(p::Point{2}, t::Triangle{2})
  a, b, c = t.vertices
  abp = signarea(a, b, p)
  bcp = signarea(b, c, p)
  cap = signarea(c, a, p)
  areas = (abp, bcp, cap)
  all(areas .≥ 0) || all(areas .≤ 0)
end
