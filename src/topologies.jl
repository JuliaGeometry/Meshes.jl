# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Topology

A data structure for constructing topological relations in a [`Mesh`](@ref).

## References

* Floriani, L. & Hui, A. 2007. [Shape representations
  based on simplicial and cell complexes]
  (https://diglib.eg.org/handle/10.2312/egst.20071055.063-087)
"""
abstract type Topology end

"""
    vertices(topology)

Return the vertices of the `topology`.
"""
vertices(t::Topology) = 1:nvertices(t)

"""
    nvertices(topology)

Return the number of vertices in the `topology`.
"""
function nvertices(::Topology) end

"""
    faces(topology, rank)

Return an iterator with `rank`-faces of the `topology`.

## Example

Consider a mesh of tetrahedrons embedded in a 3D space. We can loop over
all 3-faces (a.k.a. elements) or over all 2-faces to handle the interfaces
(i.e. triangles) between adjacent elements:

```julia
tetrahedrons = faces(topology, 3)
triangles    = faces(topology, 2)
segments     = faces(topology, 1)
```
"""
function faces(::Topology, rank) end

"""
    elements(topology)

Return the top-faces (a.k.a. elements) of the `topology`.

## Example

The elements of a volume embedded in 3D space can be tetrahedrons, hexahedrons,
or any 3-face. The elements of a surface embedded in 3D space can be triangles,
quadrangles or any 2-face.
"""
elements(t::Topology) = (element(t, i) for i in 1:nelements(t))

"""
    element(topology, ind)

Return the element of the `topology` at index `ind`.
"""
function element(::Topology, ind) end

"""
    nelements(topology)

Return the number of elements in the `topology`.
"""
function nelements(::Topology) end

"""
    facets(topology)

Return the (top-1)-faces (a.k.a. facets) of the `topology`.
"""
facets(t::Topology) = (facet(t, i) for i in 1:nfacets(t))

"""
    facet(topology, ind)

Return the facet of the `topology` at index `ind`.
"""
function facet(::Topology, ind) end

"""
    nfacets(topology)

Return the number of facets in the `topology`.
"""
function nfacets(::Topology) end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("topologies/full.jl")
include("topologies/grid.jl")
include("topologies/halfedge.jl")
