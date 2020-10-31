# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
#       v
#       ^
#       |
# 3-----------2      
# |     |     |      
# |     |     |      
# |     +---- | --> u
# |           |      
# |           |      
# 0-----------1      

"""
    Quadrangle(p1, p2, p3, p4)

A quadrangle with points `p1`, `p2`, `p3`, `p4`.
"""
struct Quadrangle{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T}
  vertices::V
end

paramdim(::Type{<:Quadrangle}) = 2

function facets(quad::Quadrangle)
  connec = connect.([(1,2),(2,3),(3,4),(4,1)], Segment)
  (materialize(c, quad.vertices) for c in connec)
end
