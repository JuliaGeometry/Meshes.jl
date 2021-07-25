# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Tetrahedron(p1, p2, p3, p4)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`.
"""
struct Tetrahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end

isconvex(::Type{<:Tetrahedron}) = true
issimplex(::Type{<:Tetrahedron}) = true

nvertices(::Type{<:Tetrahedron}) = 4
nvertices(t::Tetrahedron) = nvertices(typeof(t))

function measure(t::Tetrahedron)
  a, b, c, d = t.vertices
  abs((a - d) ⋅ ((b - d) × (c - d))) / 6
end

function edges(c::Tetrahedron)
  all_edges = ((c.vertices[1],c.vertices[2]), (c.vertices[2],c.vertices[3]),
  (c.vertices[3],c.vertices[1]), (c.vertices[1],c.vertices[4]),
  (c.vertices[2],c.vertices[4]), (c.vertices[3],c.vertices[4]))
  return (Segment([all_edges[i]...]) for i in 1:6)
end
