# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
#                    v
#                  .
#                ,/
#               /
#            2                 
#          ,/|`\               
#        ,/  |  `\             
#      ,/    '.   `\           
#    ,/       |     `\         
#  ,/         |       `\       
# 0-----------'.--------1 --> u
#  `\.         |      ,/       
#     `\.      |    ,/         
#        `\.   '. ,/           
#           `\. |/             
#              `3              
#                 `\.
#                    ` w

"""
    Tetrahedron(p1, p2, p3, p4)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`.
"""
struct Tetrahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T}
  vertices::V
end

paramdim(::Type{<:Tetrahedron}) = 3

function facets(tetra::Tetrahedron)
  connec = connect.([(1,2,3),(1,4,3),(1,4,2),(2,3,4)], Triangle)
  (materialize(c, tetra.vertices) for c in connec)
end
