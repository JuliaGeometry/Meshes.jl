# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mesh{Dim,T}

A mesh embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Mesh{Dim,T} end

"""
    faces(mesh, rank)

Return an iterator with `rank`-faces of the `mesh`.

## Example

Consider a mesh embedded in a 3D space. We can loop over all 3-faces
(a.k.a. cells) or over all 2-faces to handle the interfaces between
adjacent cells:

```julia
cells = faces(mesh, 3)
inter = faces(mesh, 2)
```
"""
function faces end

include("meshes/unstructured.jl")
