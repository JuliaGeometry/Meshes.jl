# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HalfEdge(head, elem, prev, next, half)

Stores the indices of the `head` vertex, the `prev`
and `next` edges in the left `elem`, and the opposite
`half`-edge. For some half-edges the `elem` may be
`nothing`, e.g. border edges of the mesh.

See [`HalfEdgeStructure`](@ref) for more details.
"""
mutable struct HalfEdge
  head::Int
  elem::Union{Int,Nothing}
  prev::HalfEdge
  next::HalfEdge
  half::HalfEdge
  HalfEdge(head, elem) = new(head, elem)
end

function Base.show(io::IO, e::HalfEdge)
  print(io, "HalfEdge($(e.head), $(e.elem))")
end

"""
    HalfEdgeStructure(halfedges, edgeonelem, edgeonvertex)

A data structure for orientable 2-manifolds based
on half-edges.

Two types of half-edges exist (Kettner 1999). This
implementation is the most common type that splits
the incident elements.

A vector of `halfedges` together with a vector of
`edgeonelem` and a vector of `edgeonvertex` can be
used to retrieve topolological relations in optimal
time. In this case, `edgeonvertex[i]` returns the
index of the half-edge in `halfedges` with head equal
to `i`. Similarly, `edgeonelem[i]` returns the index
of a half-edge in `halfedges` that is in the elem `i`.

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
  edgeonelem::Vector{Int}
  edgeonvertex::Vector{Int}
end

"""
    halfedges(s)

Return the half-edges of the half-edge structure `s`.
"""
halfedges(s::HalfEdgeStructure) = s.halfedges

"""
    edgeonelem(s, c)

Return a half-edge of the half-edge structure `s` on the `c`-th elem.
"""
edgeonelem(s::HalfEdgeStructure, c) = s.halfedges[s.edgeonelem[c]]

"""
    edgeonvertex(s, v)

Return the half-edge of the half-edge structure `s` for which the
head is the `v`-th index.
"""
edgeonvertex(s::HalfEdgeStructure, v) = s.halfedges[s.edgeonvertex[v]]

"""
    nelements(s)

Return the number of elements in the half-edge structure `s`.
"""
nelements(s::HalfEdgeStructure) = length(s.edgeonelem)

function boundary(c::Connectivity{<:Polygon}, ::Val{1},
                  s::HalfEdgeStructure)
  v = first(indices(c))
  e = edgeonvertex(s, v)
  n = e.next
  segments = [connect((v, n.head), Segment)]
  while n != e
    # move to next segment
    v = n.head
    n = n.next

    # add to segment list
    s = connect((v, n.head), Segment)
    push!(segments, s)
  end
  segments
end

function boundary(c::Connectivity{<:Polygon}, ::Val{0},
                  s::HalfEdgeStructure)
  v = first(indices(c))
  e = edgeonvertex(s, v)
  n = e.next
  vertices = [v]
  while n != e
    push!(vertices, n.head)
    n = n.next
  end
  vertices
end

function boundary(c::Connectivity{<:Segment}, ::Val{0},
                  s::HalfEdgeStructure)
  v = first(indices(c))
  e = edgeonvertex(s, v)
  n = e.next
  [v, n.head]
end
