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
struct Tetrahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end
