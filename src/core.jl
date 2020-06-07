
# ------------
# MESH TRAITS
# ------------
"""
    ismesh(::Type{M})
    ismesh(mesh::M)

Informs that the `mesh` implements the API.
"""
function ismesh end

"""
    ndims(::Type{M})
    ndims(mesh::M)

Number of dimensions of the spatial domain that the
`mesh` represents.
"""
function ndims end

"""
    coordtype(::Type{M})
    coordtype(mesh::M)

The machine type used to store each individual coordinate
in the embedding space of the `mesh`.
"""
function coordtype end

"""
    isstructured(::Type{M})
    isstructured(mesh::M)

Tells whether or not the `mesh` is structured, i.e. if any
element in the mesh can be located with integer indices
i, j, k, ...
"""
function isstructured end

"""
    isregular(::Type{M})
    isregular(mesh::M)

Tells whether or not the `mesh` is regular, i.e. if besides
being structured, the spacing between elements is fixed
(e.g. image).
"""
function isregular end

"""
    elements(mesh::M)

Return an iterator of `mesh` elements with known length (e.g. 1:n).
"""
function elements end

# ---------------
# ELEMENT TRAITS
# ---------------
"""
    vertices(mesh::M, elm::E)

Return the GUID of the vertices of element `elm` in the `mesh`.
"""
function vertices end

# --------------
# VERTEX TRAITS
# --------------
"""
    vcoords!(x::AbstractVector{T}, mesh::M, vert::V)

Retrieve the coordinates `x` of vertex `vert` in the `mesh`.
"""
function vcoords! end
