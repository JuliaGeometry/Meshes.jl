
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
element in the mesh can be located with integer indices.

# Example

A regular grid where each square is split into two half triangles
is structured. One can easily enumerate the elements with integer
indices and create a fast lookup scheme.
"""
function isstructured end

"""
    iscartesian(::Type{M})
    iscartesian(mesh::M)

Tells whether or not the `mesh` is Cartesian, i.e. if any
element in the mesh has "left", "right", "up" and "down"
neighbor elements. These neighbors can be located based
on the integer indices (i, j, k, ...) of the element.

# Example

Parallels and meridians in a bounding box over the Earth's
globe form a Cartesian mesh. Even though the spacing between
elements is not constant, it is still possible to navigate
the mesh with "north", "south", "east" and "west" integer
increments.
"""
function iscartesian end

"""
    isregular(::Type{M})
    isregular(mesh::M)

Tells whether or not the `mesh` is regular, i.e. if besides
being Cartesian, the spacing between elements is constant
(a.k.a. image). Regular meshes are ideal for systems with
memory constraints.

# Example

A potograph image is a regular mesh. Each pixel (or voxel in 3D)
contributes the same volume to the domain discretization.
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
