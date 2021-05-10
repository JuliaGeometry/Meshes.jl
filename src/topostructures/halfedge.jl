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

"""
    loop(halfedge)

Loop over the heads of a cycle starting at given `halfedge`.
"""
function loop(e::HalfEdge)
  n = e.next
  v = [e.head]
  while n != e
    push!(v, n.head)
    n = n.next
  end
  v
end

function Base.show(io::IO, e::HalfEdge)
  print(io, "HalfEdge($(e.head), $(e.elem))")
end

"""
    HalfEdgeStructure(elements)

A data structure for orientable 2-manifolds based
on half-edges constructed from a list of polygon
`elements`.

Two types of half-edges exist (Kettner 1999). This
implementation is the most common type that splits
the incident elements.

See also [`TopologicalStructure`](@ref).

## References

* Kettner, L. (1999). [Using generic programming for
  designing a data structure for polyhedral surfaces]
  (https://www.sciencedirect.com/science/article/pii/S0925772199000073)

### Notes

A vector of `halfedges` together with a dictionary of
`half4elem` and a dictionary of `half4vert` can be
used to retrieve topolological relations in optimal
time. In this case, `half4vert[i]` returns the index
of the half-edge in `halfedges` with head equal to `i`.
Similarly, `half4elem[i]` returns the index of a
half-edge in `halfedges` that is in the element `i`.
Additionally, a dictionary `edge4pair` returns the
index of the edge (i.e. two halves) for a given
pair of vertices.
"""
struct HalfEdgeStructure <: TopologicalStructure
  halfedges::Vector{HalfEdge}
  half4elem::Dict{Int,Int}
  half4vert::Dict{Int,Int}
  edge4pair::Dict{Tuple{Int,Int},Int}
end

function HalfEdgeStructure(elems::AbstractVector{<:Connectivity})
  @assert all(e -> paramdim(e) == 2, elems) "invalid element for half-edge structure"

  # number of vertices is the maximum index in connectivities
  nvertices = maximum(i for e in elems for i in indices(e))

  # initialization step
  half4pair = Dict{Tuple{Int,Int},HalfEdge}()
  for (e, elem) in Iterators.enumerate(elems)
    inds = collect(indices(elem))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      half4pair[(v[i], v[i+1])] = HalfEdge(v[i], e)
    end
  end

  # add missing pointers
  for elem in elems
    inds = collect(indices(elem))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      # update pointers prev and next
      he = half4pair[(v[i], v[i+1])]
      he.prev = half4pair[(v[i-1],   v[i])]
      he.next = half4pair[(v[i+1], v[i+2])]

      # if not a border element, update half
      if haskey(half4pair, (v[i+1], v[i]))
        he.half = half4pair[(v[i+1], v[i])]
      else # create half-edge for border
        be = HalfEdge(v[i+1], nothing)
        be.half = he
        he.half = be
      end
    end
  end

  # store all halfedges in a vector
  edgecount = 0
  halfedges = HalfEdge[]
  edge4pair = Dict{Tuple{Int,Int},Int}()
  for ((u, v), he) in half4pair
    if !haskey(edge4pair, (u, v))
      append!(halfedges, [he, he.half])
      edgecount += 1
      edge4pair[(u, v)] = edgecount
      edge4pair[(v, u)] = edgecount
    end
  end

  # reverse mappings
  half4elem = Dict{Int,Int}()
  half4vert = Dict{Int,Int}()
  for (e, he) in enumerate(halfedges)
    if !isnothing(he.elem) # interior half-edge
      if !haskey(half4elem, he.elem)
        half4elem[he.elem] = e
      end
      if !haskey(half4vert, he.head)
        half4vert[he.head] = e
      end
    end
  end

  HalfEdgeStructure(halfedges, half4elem, half4vert, edge4pair)
end

"""
    half4elem(e, s)

Return a half-edge of the half-edge structure `s` on the `e`-th elem.
"""
half4elem(e::Integer, s::HalfEdgeStructure) = s.halfedges[s.half4elem[e]]

"""
    half4vert(v, s)

Return the half-edge of the half-edge structure `s` for which the
head is the `v`-th index.
"""
half4vert(v::Integer, s::HalfEdgeStructure) = s.halfedges[s.half4vert[v]]

"""
    half4edge(e, s)

Return the half-edge of the half-edge structure `s` for the edge `e`.

### Notes

Always return the half-edge to the "left".
"""
half4edge(e::Integer, s::HalfEdgeStructure) = s.halfedges[2e - 1]

"""
    half4pair(uv, s)

Return the half-edge of the half-edge structure `s` for the pair of
vertices `uv`.

### Notes

Always return the half-edge to the "left".
"""
half4pair(uv::Tuple{Int,Int}, s::HalfEdgeStructure) = half4edge(s.edge4pair[uv], s)

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(s::HalfEdgeStructure) = length(s.half4vert)

function element(s::HalfEdgeStructure, ind)
  v = loop(half4elem(ind, s))
  connect(Tuple(v), Ngon{length(v)})
end

nelements(s::HalfEdgeStructure) = length(s.half4elem)

function facet(s::HalfEdgeStructure, ind)
  e = half4edge(ind, s)
  connect((e.head, e.half.head), Segment)
end

nfacets(s::HalfEdgeStructure) = length(s.halfedges) ÷ 2

# ------------
# CONVERSIONS
# ------------

Base.convert(::Type{HalfEdgeStructure}, s::TopologicalStructure) =
  HalfEdgeStructure(collect(elements(s)))
