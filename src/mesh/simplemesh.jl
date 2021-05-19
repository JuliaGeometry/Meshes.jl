# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, connec)

A simple mesh with `points` and connectivities `connec`.
The i-th face of the mesh is lazily built based on
the connectivity list `connec[i]`.
"""
struct SimpleMesh{Dim,T,TP<:Topology} <: Mesh{Dim,T}
  points::Vector{Point{Dim,T}}
  topology::TP
end

SimpleMesh(points::AbstractVector{<:Point},
           connec::AbstractVector{<:Connectivity}) =
  SimpleMesh(points, FullTopology(connec))

vertices(m::SimpleMesh) = m.points

topology(m::SimpleMesh) = m.topology
