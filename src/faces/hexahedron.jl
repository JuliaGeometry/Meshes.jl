# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
#        v
# 3----------2         
# |\     ^   |\        
# | \    |   | \       
# |  \   |   |  \      
# |   7------+---6     
# |   |  +-- |-- | -> u
# 0---+---\--1   |     
#  \  |    \  \  |     
#   \ |     \  \ |     
#    \|      w  \|     
#     4----------5     

"""
    Hexahedron(p1, p2, ..., p8)

A hexahedron with points `p1`, `p2`, ..., `p8`.
"""
struct Hexahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T}
  vertices::V
end

paramdim(::Type{<:Hexahedron}) = 3

function facets(hex::Hexahedron)
  connec = connect.([(1,2,3,4),(5,6,7,8),(1,5,6,2),
                     (2,6,7,3),(3,4,8,7),(1,5,8,4)], Quadrangle)
  (materialize(c, hex.vertices) for c in connec)
end
