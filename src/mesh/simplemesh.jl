# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, connec)

A simple mesh with `points` and connectivities `connec`.
The i-th face of the mesh is lazily built based on
the connectivity list `connec[i]`.
"""
struct SimpleMesh{Dim,T,Topology<:TopologicalStructure} <: Mesh{Dim,T}
  points::Vector{Point{Dim,T}}
  topology::Topology
end

SimpleMesh(points::AbstractVector{<:Point},
           connec::AbstractVector{<:Connectivity}) =
  SimpleMesh(points, FullStructure(connec))

==(m1::SimpleMesh, m2::SimpleMesh) =
  m1.points == m2.points && m1.topology == m2.topology

vertices(m::SimpleMesh) = m.points

faces(m::SimpleMesh, rank) =
  (materialize(c, m.points) for c in faces(m.topology, rank))

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(m::SimpleMesh, ind::Int) =
  materialize(element(m.topology, ind), m.points)

nelements(m::SimpleMesh) = nelements(m.topology)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::MIME"text/plain", m::SimpleMesh{Dim,T}) where {Dim,T}
  verts = m.points
  elms  = collect(elements(m.topology))
  nvert = length(verts)
  nelms = length(elms)
  println(io, m)
  println(io, "  $nvert vertices")
  println(io, _lines(verts, "    "))
  println(io, "  $nelms elements")
  print(  io, _lines(elms, "    "))
end

function _lines(vec, tab="  ")
  N = length(vec)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["$(tab)└─$(vec[i])" for i in 1:I]
           (N > 10 ? ["$(tab)⋮"] : [])
           ["$(tab)└─$(vec[i])" for i in J:N]]
  join(lines, "\n")
end
