# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, connec)

A simple mesh with `points` and connectivities `connec`.
The i-th face of the mesh is lazily built based on
the connectivity list `connec[i]`.

    SimpleMesh(points, topology)

Alternatively, construct a simple mesh with `points` and
a topological data structure (e.g. `HalfEdgeTopology`).

See also [`Topology`](@ref).
"""
struct SimpleMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},TP<:Topology} <: Mesh{Dim,T}
  points::V
  topology::TP
end

SimpleMesh(points::AbstractVector{<:Point},
           connec::AbstractVector{<:Connectivity}) =
  SimpleMesh(points, FullTopology(connec))

vertices(m::SimpleMesh) = m.points

topology(m::SimpleMesh) = m.topology

"""
    convert(SimpleMesh, mesh)

Convert any `mesh` to a simple mesh with explicit
list of points and [`FullTopology`](@ref).
"""
Base.convert(::Type{<:SimpleMesh}, m::Mesh) =
  topoconvert(FullTopology, m)
