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
