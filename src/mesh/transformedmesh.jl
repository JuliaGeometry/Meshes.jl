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

Base.parent(m::TransformedMesh) = m.mesh

transform(m::TransformedMesh) = m.transform

topology(m::TransformedMesh) = topology(m.mesh)

vertex(m::TransformedMesh, ind::Int) = m.transform(vertex(m.mesh, ind))

# alias to improve readability in IO methods
const TransformedGrid{Dim,T,G<:Grid{Dim,T},TR} = TransformedMesh{Dim,T,GridTopology{Dim},G,TR}

TransformedGrid(g::Grid, t::Transform) = TransformedMesh(g, t)

@propagate_inbounds Base.getindex(g::TransformedGrid{Dim}, I::CartesianIndices{Dim}) where {Dim} =
  TransformedGrid(getindex(g.mesh, I), g.transform)

function Base.summary(io::IO, g::TransformedGrid{Dim,T}) where {Dim,T}
  join(io, size(g), "×")
  print(io, " TransformedGrid{$Dim,$T}")
end
