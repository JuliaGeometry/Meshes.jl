# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HalfEdge(head, cell, prev, next, half)

Stores the indices of the `head` vertex, the `prev`
and `next` edges in the left `cell`, and the opposite
`half`-edge. For some half-edges the `cell` may be
`nothing`.

See [`HalfEdgeStructure`](@ref) for more details.
"""
mutable struct HalfEdge
  head::Int
  cell::Union{Int,Nothing}
  prev::HalfEdge
  next::HalfEdge
  half::HalfEdge
  HalfEdge(head, cell) = new(head, cell)
end

function Base.show(io::IO, e::HalfEdge)
  print(io, "HalfEdge($(e.head), $(e.cell))")
end

"""
    HalfEdgeStructure(halfedges, edgeoncell, edgeonvertex)

A data structure for orientable 2-manifolds based
on half-edges.

Two types of half-edges exist (Kettner 1999). This
implementation is the most common type that splits
the incident cells.

A vector of `halfedges` together with a vector of
`edgeoncell` and a vector of `edgeonvertex` can be
used to retrieve topolological relations in optimal
time.

Such data structure is usually constructed from another
data structure such as [`ExplicitStructure`](@ref) via
`convert` methods:

```julia
he = convert(HalfEdgeStructure, structure)
```

See also [`TopologicalStructure`](@ref).

## References

* Kettner, L. (1999). [Using generic programming for
  designing a data structure for polyhedral surfaces]
  (https://www.sciencedirect.com/science/article/pii/S0925772199000073)
"""
struct HalfEdgeStructure <: TopologicalStructure
  halfedges::Vector{HalfEdge}
  edgeoncell::Vector{Int}
  edgeonvertex::Vector{Int}
end

halfedges(s::HalfEdgeStructure) = s.halfedges
edgeoncell(s::HalfEdgeStructure, c) = s.halfedges[s.edgeoncell[c]]
edgeonvertex(s::HalfEdgeStructure, v) = s.halfedges[s.edgeonvertex[v]]
ncells(s::HalfEdgeStructure) = length(s.edgeoncell)
