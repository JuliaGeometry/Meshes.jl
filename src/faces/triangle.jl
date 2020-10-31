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

function facets(tri::Triangle)
  connec = connect.([(1,2), (2,3), (3,1)], Segment)
  (materialize(c, tri.vertices) for c in connec)
end
