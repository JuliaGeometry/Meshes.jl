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
data structure such as [`ElementListStructure`](@ref) via
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
    edgeonelem(e, s)

Return a half-edge of the half-edge structure `s` on the `e`-th elem.
"""
edgeonelem(e, s::HalfEdgeStructure) = s.halfedges[s.edgeonelem[e]]

"""
    edgeonvertex(v, s)

Return the half-edge of the half-edge structure `s` for which the
head is the `v`-th index.
"""
edgeonvertex(v, s::HalfEdgeStructure) = s.halfedges[s.edgeonvertex[v]]

# ----------------------
# TOPOLOGICAL RELATIONS
# ----------------------

function coboundary(v::Integer, ::Val{1}, s::HalfEdgeStructure)
  connect.([(v, u) for u in adjacency(v, s)], Segment)
end

function coboundary(v::Integer, ::Val{2}, s::HalfEdgeStructure)
end

function coboundary(c::Connectivity{<:Segment}, ::Val{2},
                    s::HalfEdgeStructure)
  u, v = indices(c)
  eᵤ = edgeonvertex(u, s)

  # search edge counter-clockwise
  e  = eᵤ
  while !isnothing(e.elem) && e.half.head != v
    e = e.prev.half
  end

  # search edge clockwise
  if isnothing(e.elem) && e.half.head != v
    e = eᵤ
    while !isnothing(e.elem) && e.half.head != v
      e = e.half.next
    end
  end

  # construct elements
  if isnothing(e.elem)
    [elemonedge(e.half)]
  elseif isnothing(e.half.elem)
    [elemonedge(e)]
  else
    [elemonedge(e), elemonedge(e.half)]
  end
end

function adjacency(c::Connectivity{<:Polygon}, s::HalfEdgeStructure)
end

function adjacency(c::Connectivity{<:Segment}, s::HalfEdgeStructure)
end

function adjacency(v::Integer, s::HalfEdgeStructure)
  e = edgeonvertex(v, s)
  h = e.half
  if isnothing(h.elem) # border edge
    # we are at the first arm of the star already
    # there is no need to adjust the CCW loop
  else # interior edge
    # we are at an interior edge and may need to
    # adjust the CCW loop so that it starts at
    # the first arm of the star
    n = h.next
    h = n.half
    while !isnothing(h.elem) && n != e
      n = h.next
      h = n.half
    end
    e = n
  end

  # edge e is now the first arm of the star
  # we can follow the CCW loop until we find
  # it again or hit a border edge
  p = e.prev
  n = e.next
  o = p.half
  vertices = [n.head]
  while !isnothing(o.elem) && o != e
    p = o.prev
    n = o.next
    o = p.half
    push!(vertices, n.head)
  end
  # if border edge is hit, add last arm manually
  isnothing(o.elem) && push!(vertices, o.half.head)

  vertices
end

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

element(s::HalfEdgeStructure, ind) = elemonedge(edgeonelem(ind, s))

nelements(s::HalfEdgeStructure) = length(s.edgeonelem)

# helper function
function elemonedge(e::HalfEdge)
  n = e.next
  v = [e.head]
  while n != e
    push!(v, n.head)
    n = n.next
  end
  connect(Tuple(v), Ngon{length(v)})
end
