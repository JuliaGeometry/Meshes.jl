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
    boundary(connectivity, rank, structure)

Boundary relation that maps the `connectivity` with a given
rank to a set of connectivities of smaller `rank` based on
the topological `structure`.

For example, the boundary of a `Connectivity{Triangle}` is a
set with three `Connectivity{Segment}`.
"""
boundary(c::Connectivity, rank::Integer, structure) =
  boundary(c, Val(rank), structure)

"""
    coboundary(connectivity, rank, structure)

Co-boundary relation that maps the `connectivity` with a given
rank to a set of connectivities of higher `rank` based on the
topological `structure`.

For example, the coboundary of a `Connectivity{Segment}` can
be a set with two `Connectivity{Triangle}`.
"""
coboundary(c::Connectivity, rank::Integer, structure) =
  coboundary(c, Val(rank), structure)

"""
    adjacency(connectivity, structure)

Adjacency relation that maps the `connectivity` with rank `p`
to a set of connectivities sharing a `p-1` connectivity.

For example, the adjacency of a `Connectivity{Triangle}` can
be a set of `Connectivity{Triangle}` sharing a `connectivity{Segment}`.
"""
function adjacency end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("topostructures/explicit.jl")
include("topostructures/halfedge.jl")

# conversions between structures
include("topostructures/conversions.jl")
