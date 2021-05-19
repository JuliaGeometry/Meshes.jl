# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HalfEdge(head, elem, prev, next, half)

Stores the indices of the `head` vertex, the `prev`
and `next` edges in the left `elem`, and the opposite
`half`-edge. For some half-edges the `elem` may be
`nothing`, e.g. border edges of the mesh.

See [`HalfEdgeTopology`](@ref) for more details.
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
    HalfEdgeTopology(elements)
    HalfEdgeTopology(halfedges)

A data structure for orientable 2-manifolds based on
half-edges constructed from a vector of connectivity
`elements` or from a vector of pairs of `halfedges`.

## Example

Construct half-edge topology from a list of top-faces:

```julia
elements = connect.([(1,2,3),(3,2,4,5)])
topology = HalfEdgeTopology(elements)
```

See also [`Topology`](@ref).

## References

* Kettner, L. (1999). [Using generic programming for
  designing a data structure for polyhedral surfaces]
  (https://www.sciencedirect.com/science/article/pii/S0925772199000073)

### Notes

- Two types of half-edges exist (Kettner 1999). This
  implementation is the most common type that splits
  the incident elements.

- A vector of `halfedges` together with a dictionary of
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
struct HalfEdgeTopology <: Topology
  halfedges::Vector{HalfEdge}
  half4elem::Dict{Int,Int}
  half4vert::Dict{Int,Int}
  edge4pair::Dict{Tuple{Int,Int},Int}
end

function HalfEdgeTopology(halves::AbstractVector{Tuple{HalfEdge,HalfEdge}})
  # make sure that first half-edge is in the interior
  ordered = [isnothing(h₁.elem) ? (h₂, h₁) : (h₁, h₂) for (h₁, h₂) in halves]

  # flatten pairs of half-edges into a vector
  halfedges = [half for pair in ordered for half in pair]

  # map element and vertex to a half-edge
  half4elem = Dict{Int,Int}()
  half4vert = Dict{Int,Int}()
  for (i, h) in enumerate(halfedges)
    if !isnothing(h.elem) # interior half-edge
      if !haskey(half4elem, h.elem)
        half4elem[h.elem] = i
      end
      if !haskey(half4vert, h.head)
        half4vert[h.head] = i
      end
    end
  end

  # map pair of vertices to an edge (i.e. two halves)
  edge4pair = Dict{Tuple{Int,Int},Int}()
  for (i, (h₁, h₂)) in enumerate(ordered)
    u, v = h₁.head, h₂.head
    edge4pair[(u, v)] = i
    edge4pair[(v, u)] = i
  end

  HalfEdgeTopology(halfedges, half4elem, half4vert, edge4pair)
end

function HalfEdgeTopology(elems::AbstractVector{<:Connectivity})
  @assert all(e -> paramdim(e) == 2, elems) "invalid element for half-edge topology"

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

  # save halfedges in a vector of pairs
  halves  = Vector{Tuple{HalfEdge,HalfEdge}}()
  visited = Set{Tuple{Int,Int}}()
  for ((u, v), he) in half4pair
    if (u, v) ∉ visited
      push!(halves,  (he, he.half))
      push!(visited, (u, v))
      push!(visited, (v, u))
    end
  end

  HalfEdgeTopology(halves)
end

"""
    half4elem(e, t)

Return a half-edge of the half-edge topology `t` on the `e`-th elem.
"""
half4elem(e::Integer, t::HalfEdgeTopology) = t.halfedges[t.half4elem[e]]

"""
    half4vert(v, t)

Return the half-edge of the half-edge topology `t` for which the
head is the `v`-th index.
"""
half4vert(v::Integer, t::HalfEdgeTopology) = t.halfedges[t.half4vert[v]]

"""
    half4edge(e, t)

Return the half-edge of the half-edge topology `t` for the edge `e`.

### Notes

Always return the half-edge to the "left".
"""
half4edge(e::Integer, t::HalfEdgeTopology) = t.halfedges[2e - 1]

"""
    half4pair(uv, t)

Return the half-edge of the half-edge topology `t` for the pair of
vertices `uv`.

### Notes

Always return the half-edge to the "left".
"""
half4pair(uv::Tuple{Int,Int}, t::HalfEdgeTopology) = half4edge(t.edge4pair[uv], t)

"""
    edge4pair(uv, t)

Return the edge of the half-edge topology `t` for the pair of vertices `uv`.
"""
edge4pair(uv, t) = t.edge4pair[uv]

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::HalfEdgeTopology) = length(t.half4vert)

function faces(t::HalfEdgeTopology, rank)
  if rank == 2
    elements(t)
  elseif rank == 1
    facets(t)
  else
    throw(ArgumentError("invalid rank for half-edge topology"))
  end
end

function element(t::HalfEdgeTopology, ind)
  v = loop(half4elem(ind, t))
  connect(Tuple(v))
end

nelements(t::HalfEdgeTopology) = length(t.half4elem)

function facet(t::HalfEdgeTopology, ind)
  e = half4edge(ind, t)
  connect((e.head, e.half.head))
end

nfacets(t::HalfEdgeTopology) = length(t.halfedges) ÷ 2

# ------------
# CONVERSIONS
# ------------

Base.convert(::Type{<:HalfEdgeTopology}, t::Topology) =
  HalfEdgeTopology(collect(elements(t)))
