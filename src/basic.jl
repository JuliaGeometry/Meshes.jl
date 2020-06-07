"""
    nelms(mesh::M)

Number of elements in the `mesh`.
"""
nelms(mesh::M) where M = length(elements(mesh))

"""
   nverts(mesh::M, elm::E)

Number of vertices in element `elm` of the `mesh`.
"""
nverts(mesh::M, elm::E) where {M,E} = length(vertices(mesh, elm))

"""
    coordbuff(mesh)

A buffer for storing coordinates in a reference system.
"""
coordbuff(mesh::Type{M}) where M = MVector{ndims(mesh),coordtype(mesh)}(undef)

"""
    vcoords!(X, mesh::M, verts::AbstractVector{V})

Retrieve coordinates `X` of vertices `verts` in `mesh`
as columns of the matrix.
"""
function vcoords!(X::AbstractMatrix, mesh::M,
                  verts::AbstractVector{V}) where {M,V}
  for j in 1:length(verts)
    @inbounds vcoords!(view(X,:,j), mesh, verts[j])
  end
end

"""
    vcoords(mesh::M, vert::V)

Return coordinates of vertex `vert` in `mesh`.
"""
function vcoords(mesh::M, vert::V) where {M,V}
  x = coordbuff(M)
  vcoords!(x, mesh, vert)
  x
end

"""
    vcoords(mesh::M, verts::AbstractVector{V})

Return coordinates of vertices `verts` in `mesh`.
"""
function vcoords(mesh::M, verts::AbstractVector{V}) where {M,V}
  X = Matrix{coordtype(M)}(undef, ndims(M), length(verts))
  vcoords!(X, mesh, verts)
  X
end
