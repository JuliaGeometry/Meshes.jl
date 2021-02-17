# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

"""
    signarea(triangle)

Compute the signed area of `triangle`.
"""
function signarea(t::Triangle{2})
  vs = vertices(t)
  signarea(vs[1], vs[2], vs[3])
end

"""
    sideof(point, segment)

Determines on which side of the oriented `segment`
the `point` lies. Possible results are `:LEFT`,
`:RIGHT` or `:ON` the segment.
"""
function sideof(p::Point{2,T}, s::Segment{2,T}) where {T}
  a, b = vertices(s)
  area = signarea(p, a, b)
  ifelse(area > atol(T), :LEFT, ifelse(area < -atol(T), :RIGHT, :ON))
end
