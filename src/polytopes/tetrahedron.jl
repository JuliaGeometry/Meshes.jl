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
  A, B, C, D = t.vertices
  abs((A - D) ⋅ ((B - D) × (C - D))) / 6
end

function boundary(t::Tetrahedron)
  indices = [(1,2,3),(2,1,4),(2,4,3),(3,4,1)]
  SimpleMesh(t.vertices, connect.(indices))
end