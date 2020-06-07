# -----------------
# GEOMETRIC TRAITS
# -----------------
"""
    ecoords!(x::AbstractVector{T}, mesh::M, elm::M)

Retrieve the coordinates `x` of the centroid of the element
`elm` in `mesh`. Default to mean of vertices of `elm`.
"""
function ecoords!(x::AbstractVector, mesh::M, elm::E) where {M,E}
  X = vcoords(mesh, vertices(mesh, elm))
  nd, nv = size(X)
  for i in 1:nd
    @inbounds x[i] = sum(view(X,i,:)) / nv
  end
end

"""
    ecoords!(X, mesh, elms)

Retrieve the coordinates of the centroids of elements
`elms` in the `mesh` as columns of `X`.
"""
function ecoords!(X::AbstractMatrix, mesh::M,
                  elms::AbstractVector{E}) where {M,E}
  for j in 1:length(elms)
    ecoords!(view(X,:,j), mesh, elms[j])
  end
end

"""
    ecoords(mesh, elm)

Return the coordinates of of the centroid of the element
`elm` in `mesh`.
"""
function ecoords(mesh::M, elm::E) where {M,E}
  x = coordbuff(M)
  ecoords!(x, mesh, elm)
  x
end

"""
    ecoords(mesh, elms)

Return the coordinates of the centroids of elements
`elms` in `mesh`.
"""
function ecoords(mesh::M, elms::AbstractVector{E}) where {M,E}
  X = Matrix{coordtype(M)}(undef, ndims(M), length(elms))
  ecoords!(X, mesh, elms)
  X
end

"""
    ecoords(mesh)

Return the coordinates of the centroids of all elements in `mesh`.
"""
ecoords(mesh::M) where M = ecoords(mesh, elements(mesh))

"""
    volume(mesh::M, elm::E)

Volume of `mesh` element `elm`.
"""
function volume end

"""
    volume(mesh::M, elms::AbstractVector{E})

Sum of volumes for various elements `elms` in `mesh`.
"""
volume(mesh::M, elms::AbstractVector{E}) where {M,E} =
  sum(volume(mesh, elm) for elm in elms)

"""
    volume(mesh::M)

Sum of volumes for all elements in `mesh`.
"""
volume(mesh::M) where M = volume(mesh, elements(mesh))
