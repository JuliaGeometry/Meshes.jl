# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TopologicalRelation

A topological relation between faces of a [`Mesh`](@ref) implemented
for a given [`Topology`](@ref).

## References

* Floriani, L. & Hui, A. 2007. [Shape representations
  based on simplicial and cell complexes]
  (https://diglib.eg.org/handle/10.2312/egst.20071055.063-087)
"""
abstract type TopologicalRelation end

"""
    relation(face)

Evaluate the topological `relation` at a `face` represented by an index.
"""
function (::TopologicalRelation)(::Integer) end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("toporelations/boundary.jl")
include("toporelations/coboundary.jl")
include("toporelations/adjacency.jl")
