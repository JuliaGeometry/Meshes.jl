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

topology(m::SimpleMesh) = m.topology
