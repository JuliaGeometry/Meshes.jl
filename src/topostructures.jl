# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TopologicalStructure

A data structure for constructing topological relations
in a [`Mesh`](@ref).

## References

* Floriani, L. & Hui, A. 2007. [Shape representations
  based on simplicial and cell complexes]
  (https://diglib.eg.org/handle/10.2312/egst.20071055.063-087)
"""
abstract type TopologicalStructure end

"""
    vertices(structure)

Return the vertices of the topological `structure`.
"""
vertices(s::TopologicalStructure) = 1:nvertices(s)

"""
    nvertices(structure)

Return the number of vertices in topological `structure`.
"""
function nvertices(::TopologicalStructure) end

"""
    faces(structure, rank)

Return an iterator with `rank`-faces of the topological `structure`.

## Example

Consider a mesh of tetrahedrons embedded in a 3D space. We can loop over
all 3-faces (a.k.a. elements) or over all 2-faces to handle the interfaces
(i.e. triangles) between adjacent elements:

```julia
tetrahedrons = faces(structure, 3)
triangles    = faces(structure, 2)
segments     = faces(structure, 1)
```
"""
function faces(::TopologicalStructure, rank) end

"""
    elements(structure)

Return the top-faces (a.k.a. elements) of the topological `structure`.

## Example

The elements of a volume embedded in 3D space can be tetrahedrons, hexahedrons,
or any 3-face. The elements of a surface embedded in 3D space can be triangles,
quadrangles or any 2-face.
"""
elements(structure::TopologicalStructure) =
  (element(structure, ind) for ind in 1:nelements(structure))

"""
    element(structure, ind)

Return the element of the topological `structure` at index `ind`.
"""
function element(::TopologicalStructure, ind) end

"""
    nelements(structure)

Return the number of elements in the topological `structure`.
"""
function nelements(::TopologicalStructure) end

"""
    facets(structure)

Return the (top-1)-faces (a.k.a. facets) of the topological `structure`.
"""
facets(structure::TopologicalStructure) =
  (facet(structure, ind) for ind in 1:nfacets(structure))

"""
    facet(structure, ind)

Return the facet of the topological `structure` at index `ind`.
"""
function facet(::TopologicalStructure, ind) end

"""
    nfacets(structure)

Return the number of facets in the topological `structure`.
"""
function nfacets(::TopologicalStructure) end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("topostructures/full.jl")
include("topostructures/halfedge.jl")
