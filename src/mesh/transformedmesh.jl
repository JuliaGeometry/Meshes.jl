# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TransformedMesh(mesh, transform)

Lazy representation of a geometric `transform` applied to a `mesh`.
"""
struct TransformedMesh{Dim,T,TP<:Topology,M<:Mesh{Dim,T,TP},TR<:Transform} <: Mesh{Dim,T,TP}
  mesh::M
  transform::TR
end

# specialize constructor to avoid deep structures
TransformedMesh(m::TransformedMesh, t::Transform) = TransformedMesh(m.mesh, m.transform → t)

# alias to improve readability in IO methods
const TransformedGrid{Dim,T} = TransformedMesh{Dim,T,GridTopology{Dim}}

TransformedGrid(g::Grid, t::Transform) = TransformedMesh(g, t)

topology(m::TransformedMesh) = topology(m.mesh)

vertex(m::TransformedMesh, ind::Int) = m.transform(vertex(m.mesh, ind))

Base.getindex(g::TransformedGrid{Dim}, I::CartesianIndices{Dim}) where {Dim} =
  TransformedGrid(getindex(g.mesh, I), g.transform)
