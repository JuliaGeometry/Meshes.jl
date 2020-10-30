# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
# v
# ^                 
# |                 
# 2                 
# |`\               
# |  `\             
# |    `\           
# |      `\         
# |        `\       
# 0----------1 --> u

"""
    Triangle(p1, p2, p3)

A triangle with points `p1`, `p2`, `p3`.
"""
struct Triangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T,2}
  vertices::V
end

function volume(tri::Triangle)
  A, B, C = tri.vertices
  abs((B - A) × (C - A)) / 2
end

function Base.in(P::Point, tri::Triangle)
  A, B, C = tri.vertices
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

function facets(tri::Triangle)
  connec = connect.([(1,2), (2,3), (3,1)], Segment)
  (materialize(c, tri.vertices) for c in connec)
end
