# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GMSH

Ordering convention of the GMSH project:
https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
"""
struct GMSH <: OrderingConvention end

# Triangle
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

connectivities(::Type{<:Triangle}, ::Val{1}, ::Type{GMSH}) =
  connect.([(1, 2), (2, 3), (3, 1)], Segment)

# Quadrangle
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

connectivities(::Type{<:Quadrangle}, ::Val{1}, ::Type{GMSH}) =
  connect.([(1, 2), (2, 3), (3, 4), (4, 1)], Segment)

# Pyramid
#                4
#              ,/|\
#            ,/ .'|\
#          ,/   | | \
#        ,/    .' | `.
#      ,/      |  '.  \
#    ,/       .' w |   \
#  ,/         |  ^ |    \
# 0----------.'--|-3    `.
#  `\        |   |  `\    \
#    `\     .'   +----`\ - \ -> v
#      `\   |    `\     `\  \
#        `\.'      `\     `\`
#           1----------------2
#                     `\
#                        u

function connectivities(::Type{<:Pyramid}, ::Val{2}, ::Type{GMSH})
  base = connect((1,2,3,4), Quadrangle)
  side = connect.([(1,2,5),(2,3,5),(3,4,5),(4,1,5)], Triangle)
  [base; side]
end

# Tetrahedron
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

connectivities(::Type{<:Tetrahedron}, ::Val{2}, ::Type{GMSH}) =
  connect.([(1,2,3),(1,4,3),(1,4,2),(2,3,4)], Triangle)

# Hexahedron
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

connectivities(::Type{<:Hexahedron}, ::Val{2}, ::Type{GMSH}) =
  connect.([(1,2,3,4),(5,6,7,8),(1,5,6,2),
            (2,6,7,3),(3,4,8,7),(1,5,8,4)], Quadrangle)
