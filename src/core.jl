
#--------------
# MESH TRAITS
#--------------
"""
    ndims(mesh)

Number of dimensions of the `mesh`.
"""
ndims(::Type{M}) where M = @error "not implemented"

"""
    ctype(mesh)

Coordinate type of the `mesh`, i.e coordinates are
stored as static vectors with entries of this type.
"""
ctype(::Type{M}) where M = @error "not implemented"

"""
    isstructured(mesh)

Tells whether or not the `mesh` is structured, i.e. if any element
in the `mesh` can be located with integer indices i, j, k, ...
"""
isstructured(::Type{M}) where M = Val{false}()

"""
    isregular(mesh)

Tells whether or not the `mesh` is regular, i.e. if besides being
structured, the spacing between elements is fixed (e.g. image).

### Notes

Regular meshes have low memory storage in Cartesian spaces.
"""
isregular(::Type{M}) where M = Val{false}()

"""
    elements(mesh)

An iterator of `mesh` elements with known length (e.g. 1:n).
"""
elements(mesh::M) where M = @error "not implemented"

#----------------------
# MESH ELEMENT TRAITS
#----------------------
"""
    vertices(mesh, elm)

Return the vertices of element `elm` in the `mesh`.
"""
vertices(mesh::M, elm::E) where {M,E} = @error "not implemented"

#------------------------
# ELEMENT VERTEX TRAITS
#------------------------
"""
    coords!(x, mesh, elm, vert)

Set the coordinates `x` of vertex `vert` in element `elm` in the `mesh`.
"""
coords!(x::AbstractVector, mesh::M, elm::E, vert::V) where {M,E,V} =
  @error "not implemented"
