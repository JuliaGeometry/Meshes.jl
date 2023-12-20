# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TransformedMesh(mesh, transform)

Lazy representation of a geometric `transform` applied to a `mesh`.
"""
struct TransformedMesh{Dim,T,M<:Mesh{Dim,T},TR<:Transform} <: Mesh{Dim,T}
  mesh::M
  transform::TR
end

TransformedMesh(mesh::TransformedMesh, transform::Transform) = TransformedMesh(mesh.mesh, mesh.transform → transform)

topology(m::TransformedMesh) = topology(m.mesh)

vertex(m::TransformedMesh, ind::Int) = m.transform(vertex(m.mesh, ind))
