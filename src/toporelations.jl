# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TopologicalRelation

A topological relation between faces of a [`Mesh`](@ref) implemented
for a given [`Topology`](@ref).

An object implementing this trait is a functor that can be evaluated
at an integer index representing the face.

## Examples

```julia
# create boundary relation mapping
# 2-faces to 0-faces (i.e. vertices)
∂ = Boundary{2,0}(topology)

# list of vertices for first face
∂(1)
```

## References

* Floriani, L. & Hui, A. 2007. [Shape representations
  based on simplicial and cell complexes]
  (https://diglib.eg.org/handle/10.2312/egst.20071055.063-087)
"""
abstract type TopologicalRelation end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("toporelations/boundary.jl")
include("toporelations/coboundary.jl")
include("toporelations/adjacency.jl")
