# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TransformedMesh(mesh, transform)

Lazy representation of a coordinate `transform` applied to a `mesh`.
"""
struct TransformedMesh{M<:Manifold,C<:CRS,TP<:Topology,MS<:Mesh,TR<:Transform} <: Mesh{M,C,TP}
  mesh::MS
  transform::TR

  function TransformedMesh{M,C}(mesh::MS, transform::TR) where {M<:Manifold,C<:CRS,MS<:Mesh,TR<:Transform}
    TP = typeof(topology(mesh))
    new{M,C,TP,MS,TR}(mesh, transform)
  end
end

function TransformedMesh(m::Mesh, t::Transform)
  p = t(vertex(m, 1))
  TransformedMesh{manifold(p),crs(p)}(m, t)
end

# specialize constructor to avoid deep structures
TransformedMesh(m::TransformedMesh, t::Transform) = TransformedMesh(m.mesh, m.transform → t)

Base.parent(m::TransformedMesh) = m.mesh

transform(m::TransformedMesh) = m.transform

topology(m::TransformedMesh) = topology(m.mesh)

vertex(m::TransformedMesh, ind::Int) = m.transform(vertex(m.mesh, ind))

==(m₁::TransformedMesh, m₂::TransformedMesh) = m₁.transform == m₂.transform && m₁.mesh == m₂.mesh

# alias to improve readability in IO methods
const TransformedGrid{M<:Manifold,C<:CRS,Dim,G<:Grid,TR<:Transform} = TransformedMesh{M,C,GridTopology{Dim},G,TR}

TransformedGrid(g::Grid, t::Transform) = TransformedMesh(g, t)

@propagate_inbounds Base.getindex(g::TransformedGrid, I::CartesianIndices) =
  TransformedGrid(getindex(g.mesh, I), g.transform)

function Base.summary(io::IO, g::TransformedGrid)
  join(io, size(g), "×")
  print(io, " TransformedGrid")
end
