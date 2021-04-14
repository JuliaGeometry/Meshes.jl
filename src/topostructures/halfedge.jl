# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HalfEdge(head, cell, prev, next, half)

Stores the indices of the `head` vertex, the `prev`
and `next` edges in the left `cell`, and the opposite
`half`-edge. For some half-edges the `cell` may be
`nothing`, e.g. border edges of the mesh.

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
time. In this case, `edgeonvertex[i]` returns the
index of the half-edge in `halfedges` with head equal
to `i`. Similarly, `edgeoncell[i]` returns the index
of a half-edge in `halfedges` that is in the cell `i`.

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

"""
    halfedges(s)

Return the half-edges of the half-edge structure `s`.
"""
halfedges(s::HalfEdgeStructure) = s.halfedges

"""
    edgeoncell(s, c)

Return a half-edge of the half-edge structure `s` on the `c`-th cell.
"""
edgeoncell(s::HalfEdgeStructure, c) = s.halfedges[s.edgeoncell[c]]

"""
    edgeonvertex(s, v)

Return the half-edge of the half-edge structure `s` for which the
head is the `v`-th index.
"""
edgeonvertex(s::HalfEdgeStructure, v) = s.halfedges[s.edgeonvertex[v]]

"""
    ncells(s)

Return the number of cells in the half-edge structure `s`.
"""
ncells(s::HalfEdgeStructure) = length(s.edgeoncell)
