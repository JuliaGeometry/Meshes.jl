"""
    boundary(mesh)

An iterator for the elements in boundary of the `mesh`.
"""
boundary(mesh::M) where M = @error "not implemented"

interior(mesh::M) where M = setdiff(elements(mesh), boundary(mesh))

"""
    center!(x, mesh, elm)

Set the coordinates `x` of the centroid of the `mesh` `elm` in place.
"""
center!(x::AbstractVector, mesh::M, elm::E) where {M,E} =
  @error "not implemented"

"""
    center!(X, mesh, elms)

Set the coordinates of the centroids of elements `elms` in
the `mesh` as columns of `X`.
"""
function center!(X::AbstractMatrix, mesh::M,
                 elms::AbstractVector) where M
  for j in 1:length(elms)
    center!(view(X,:,j), mesh, elms[j])
  end
end

"""
    center(mesh, elm)

Allocating version of `center!`.
"""
function center(mesh::M, elm::E) where {M,E}
  x = cbuff(M)
  center!(x, mesh, elm)
  x
end

"""
    center(mesh, elms)

Allocating version of `center!`.
"""
function center(mesh::M, elms::AbstractVector) where M
  X = Matrix{ctype(M)}(undef, ndims(M), length(elms))
  center!(X, mesh, elms)
  X
end

volume(mesh::M, elm::E) where {M,E} = @error "not implemented"

volume(mesh::M, elms::AbstractVector) where M =
  sum(volume(mesh, elm) for elm in elms)

volume(mesh::M) where M = volume(mesh, elements(mesh))
