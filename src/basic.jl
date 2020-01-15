"""
    nelms(mesh)

Number of elements in the `mesh`.
"""
nelms(mesh::M) where M = length(elements(mesh))

"""
   nverts(mesh, elm)

Number of vertices in element `elm` of the `mesh`.
"""
nverts(mesh::M, elm::E) where {M,E} = length(vertices(mesh, elm))

"""
    cbuff(mesh)

A buffer for storing coordinates in a reference system.
"""
cbuff(m::Type{M}) where M = MVector{ndims(m),ctype(m)}(undef)

"""
    coords!(X, mesh, elm, verts)

Set coordinates `X` of vertices `verts` of element `elm`
in `mesh` as columns of the matrix.
"""
function coords!(X::AbstractMatrix, mesh::M, elm::E,
                 verts::AbstractVector) where {M,E}
  for j in 1:length(verts)
    coords!(view(X,:,j), mesh, elm, verts[j])
  end
end

"""
    coords(mesh, elm, vert)

Allocating version of `coords!`.
"""
function coords(mesh::M, elm::E, vert::V) where {M,E,V}
  x = cbuff(M)
  coords!(x, mesh, elm, vert)
  x
end

"""
    coords(mesh, elm, verts)

Allocating version of `coords!`.
"""
function coords(mesh::M, elm::E, verts::AbstractVector) where {M,E}
  X = Matrix{ctype(M)}(undef, ndims(M), nelms(mesh))
  coords!(X, mesh, elm, verts)
  X
end
