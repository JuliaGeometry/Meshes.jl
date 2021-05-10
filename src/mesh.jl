# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mesh{Dim,T}

A mesh embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Mesh{Dim,T} <: Domain{Dim,T} end

"""
    vertices(mesh)

Return the vertices of the `mesh`.
"""
function vertices end

"""
    topology(mesh)

Return the topological structure of the `mesh`.
"""
function topology end

"""
    faces(mesh, rank)

Return an iterator with `rank`-faces of the `mesh`.

    faces(mesh)

Return an iterator with all faces of the `mesh`.

## Example

Consider a mesh of tetrahedrons embedded in a 3D space. We can loop over
all 3-faces (a.k.a. cells) or over all 2-faces to handle the interfaces
between adjacent cells:

```julia
tetrahedrons = faces(mesh, 3)
triangles = faces(mesh, 2)
segments = faces(mesh, 1)
```
"""
faces(m::Mesh{Dim}) where {Dim} = Iterators.flatten(faces(m, r) for r in 1:Dim)

"""
    elements(mesh)

Return the faces of the mesh with rank equal to the embedding dimension.

## Example

A 2D surface embedded in 3D space has no elements. When embedded in 2D space
the elements can be triangles, quadrangles or any 2-face.

The elements of a volume embedded in 3D space can be tetrahedrons, hexahedrons,
of any 3-face.
"""
elements(m::Mesh{Dim}) where {Dim} = faces(m, Dim)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("mesh/cartesiangrid.jl")
include("mesh/simplemesh.jl")
