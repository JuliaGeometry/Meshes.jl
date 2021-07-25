# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Hexahedron(p1, p2, ..., p8)

A hexahedron with points `p1`, `p2`, ..., `p8`.
"""
struct Hexahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end

nvertices(::Type{<:Hexahedron}) = 8
nvertices(h::Hexahedron) = nvertices(typeof(h))
function edges(c::Hexahedron)
  all_edges = ((c.vertices[1],c.vertices[2]), (c.vertices[2],c.vertices[3]),
  (c.vertices[3],c.vertices[4]), (c.vertices[4],c.vertices[1]),
  (c.vertices[1],c.vertices[5]), (c.vertices[2],c.vertices[6]),
  (c.vertices[3],c.vertices[7]), (c.vertices[4],c.vertices[8]),
  (c.vertices[5],c.vertices[6]), (c.vertices[6],c.vertices[7]),
  (c.vertices[7],c.vertices[8]), (c.vertices[8],c.vertices[5]))
  (Segment([all_edges[i]...]) for i in 1:12)
end
