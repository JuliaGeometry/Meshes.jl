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
    vertcoords!(X, mesh::M, verts::AbstractVector{V})

Retrieve coordinates `X` of vertices `verts` in `mesh`
as columns of the matrix.
"""
function vertcoords!(X::AbstractMatrix, mesh::M,
                  verts::AbstractVector{V}) where {M,V}
  for j in 1:length(verts)
    @inbounds vertcoords!(view(X,:,j), mesh, verts[j])
  end
end

"""
    vertcoords(mesh::M, vert::V)

Return coordinates of vertex `vert` in `mesh`.
"""
function vertcoords(mesh::M, vert::V) where {M,V}
  x = coordbuff(M)
  vertcoords!(x, mesh, vert)
  x
end

"""
    vertcoords(mesh::M, verts::AbstractVector{V})

Return coordinates of vertices `verts` in `mesh`.
"""
function vertcoords(mesh::M, verts::AbstractVector{V}) where {M,V}
  X = Matrix{coordtype(M)}(undef, ndims(M), length(verts))
  vertcoords!(X, mesh, verts)
  X
end
