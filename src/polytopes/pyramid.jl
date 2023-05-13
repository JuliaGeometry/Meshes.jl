# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
struct Pyramid{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end

nvertices(::Type{<:Pyramid}) = 5
nvertices(p::Pyramid) = nvertices(typeof(p))

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(p.vertices, connect.(indices))
end
