module Meshes

using StaticArrays: MVector

#--------------
# MESH TRAITS
#--------------
ndims(::Type{M}) where M = @error "not implemented"

ctype(::Type{M}) where M = @error "not implemented"

cbuff(m::Type{M}) where M = MVector{ndims(m),ctype(m)}(undef)

isstructured(::Type{M}) where M = Val{false}()

isregular(::Type{M}) where M = Val{false}()

COMPILE_TIME_TRAITS = [:ndims, :ctype, :cbuff, :isstructured, :isregular]

# default versions for mesh instances
for TRAIT in COMPILE_TIME_TRAITS
  @eval $TRAIT(::M) where M = $TRAIT(M)
end

elements(mesh::M) where M = @error "not implemented"

boundary(mesh::M) where M = @error "not implemented"

interior(mesh::M) where M = setdiff(elements(mesh), boundary(mesh))

nelms(mesh::M) where M = length(elements(mesh))

#----------------------
# MESH ELEMENT TRAITS
#----------------------
coords!(x::AbstractVector, mesh::M, elm::E) where {M,E} =
  @error "not implemented"

function coords!(X::AbstractMatrix, mesh::M,
                 elms::AbstractVector) where M
  for j in 1:length(elms)
    coords!(view(X,:,j), mesh, elms[j])
  end
end

function coords(mesh::M, elm::E) where {M,E}
  x = cbuff(M)
  coords!(x, mesh, elm)
  x
end

function coords(mesh::M, elms::AbstractVector) where M
  X = Matrix{ctype(M)}(undef, ndims(M), length(elms))
  coords!(X, mesh, elms)
  X
end

coords(mesh::M) where M = coords(mesh, elements(mesh))

volume(mesh::M, elm::E) where {M,E} = @error "not implemented"

volume(mesh::M, elms::AbstractVector) where M =
  sum(volume(mesh, elm) for elm in elms)

volume(mesh::M) where M = volume(mesh, elements(mesh))

vertices(mesh::M, elm::E) where {M,E} = @error "not implemented"

#------------------------
# ELEMENT VERTEX TRAITS
#------------------------
coords!(x::AbstractVector, mesh::M, elm::E, vert::V) where {M,E,V} =
  @error "not implemented"

function coords!(X::AbstractMatrix, mesh::M, elm::E,
                 verts::AbstractVector) where {M,E}
  for j in 1:length(verts)
    coords!(view(X,:,j), mesh, elm, verts[j])
  end
end

function coords(mesh::M, elm::E, vert::V) where {M,E,V}
  x = cbuff(M)
  coords!(x, mesh, elm, vert)
  x
end

function coords(mesh::M, elm::E, verts::AbstractVector) where {M,E}
  X = Matrix{ctype(M)}(undef, ndims(M), nelms(mesh))
  coords!(X, mesh, elm, verts)
  X
end

end # module
