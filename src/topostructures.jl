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

# ----------------------
# TOPOLOGICAL RELATIONS
# ----------------------

"""
    boundary(connectivity, rank, structure)

Boundary relation that maps the `connectivity` with a given
rank to a set of connectivities of smaller `rank` based on
the topological `structure`.

For example, the boundary of a `Connectivity{Triangle}` is a
set with three `Connectivity{Segment}`.
"""
boundary(connectivity, rank, structure::TopologicalStructure) =
  boundary(connectivity, Val(rank), structure)

function boundary(c::Connectivity{<:Polygon}, ::Val{1}, ::TopologicalStructure)
  v = CircularVector(collect(indices(c)))
  [connect((v[i], v[i+1]), Segment) for i in 1:length(v)]
end

function boundary(c::Connectivity{<:Polytope}, ::Val{0}, ::TopologicalStructure)
  collect(indices(c))
end

"""
    coboundary(connectivity, rank, structure)

Co-boundary relation that maps the `connectivity` with a given
rank to a set of connectivities of higher `rank` based on the
topological `structure`.

For example, the coboundary of a `Connectivity{Segment}` can
be a set with two `Connectivity{Triangle}`.
"""
coboundary(connectivity, rank, structure::TopologicalStructure) =
  coboundary(connectivity, Val(rank), structure)

"""
    adjacency(connectivity, structure)

Adjacency relation that maps the `connectivity` with rank `p`
to a set of connectivities sharing a `p-1` connectivity.

For example, the adjacency of a `Connectivity{Triangle}` can
be a set of `Connectivity{Triangle}` sharing a `Connectivity{Segment}`.
"""
function adjacency end

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

"""
    faces(structure, rank)

Return an iterator with `rank`-faces of the topological `structure`.

## Example

Consider a mesh of tetrahedrons embedded in a 3D space. We can loop over
all 3-faces (a.k.a. elements) or over all 2-faces to handle the interfaces
(i.e. triangles) between adjacent elements:

```julia
tetrahedrons = faces(structure, 3)
triangles = faces(structure, 2)
segments = faces(structure, 1)
```
"""
faces(structure::TopologicalStructure, rank) = faces(structure, Val(rank))

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

# ----------------
# IMPLEMENTATIONS
# ----------------

include("topostructures/full.jl")
include("topostructures/halfedge.jl")
