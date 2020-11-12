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
struct Hexahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end
