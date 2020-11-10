# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
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

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
struct Pyramid{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end

function facets(pyr::Pyramid)
  base = connect((1,2,3,4), Quadrangle)
  side = connect.([(1,2,5),(2,3,5),(3,4,5),(4,1,5)], Triangle)
  connec = [base; side]
  (materialize(c, pyr.vertices) for c in connec)
end
