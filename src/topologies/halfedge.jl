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
    HalfEdgeTopology(elements; sort=true)
    HalfEdgeTopology(halfedges)

A data structure for orientable 2-manifolds based on
half-edges constructed from a vector of connectivity
`elements` or from a vector of pairs of `halfedges`.

The option `sort` can be used to sort the elements in
adjacent-first order in case of inconsistent orientation
(i.e. mix of clockwise and counter-clockwise).

## Examples

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

Two types of half-edges exist (Kettner 1999). This
implementation is the most common type that splits
the incident elements.

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

If the `elements` of the mesh already have consistent
orientation, then the `sort` option can be disabled
for maximum performance.
"""
struct HalfEdgeTopology <: Topology
  halfedges::Vector{HalfEdge}
  half4elem::Dict{Int,Int}
  half4vert::Dict{Int,Int}
  edge4pair::Dict{Tuple{Int,Int},Int}
end

function HalfEdgeTopology(halves::AbstractVector{Tuple{HalfEdge,HalfEdge}}, nelems::Int)
  halfedges = Vector{HalfEdge}(undef, 2 * length(halves))
  edge4pair = Dict{Tuple{Int,Int},Int}()
  half4elem = Dict{Int,Int}()
  half4vert = Dict{Int,Int}()
  sizehint!(edge4pair, length(halves))
  sizehint!(half4elem, nelems)

  # flatten pairs of half-edges into a vector
  for (i, (h₁, h₂)) in enumerate(halves)
    # make sure that first half-edge is in the interior
    (h₁, h₂) = isnothing(h₁.elem) ? (h₂, h₁) : (h₁, h₂)

    j = 2i - 1
    halfedges[j] = h₁
    halfedges[j + 1] = h₂

    # map element and vertex to a half-edge
    for (k, h) in zip(j:(j + 1), (h₁, h₂))
      if !isnothing(h.elem) # interior half-edge
        get!(half4elem, h.elem, k)
        get!(half4vert, h.head, k)
      end
    end

    # map pair of vertices to an edge (i.e. two halves)
    u, v = h₁.head, h₂.head
    edge4pair[minmax(u, v)] = i
  end

  HalfEdgeTopology(halfedges, half4elem, half4vert, edge4pair)
end

function HalfEdgeTopology(elems::AbstractVector{<:Connectivity}; sort=true)
  assertion(all(e -> paramdim(e) == 2, elems), "invalid element for half-edge topology")

  # sort elements to make sure that they
  # are traversed in adjacent-first order
  adjelems = sort ? adjsort(elems) : elems
  eleminds = sort ? indexin(adjelems, elems) : 1:length(elems)

  # start assuming that all elements are
  # oriented consistently as CCW
  CCW = trues(length(adjelems))

  # initialize with first element
  half4pair = Dict{Tuple{Int,Int},HalfEdge}()
  elem = first(adjelems)
  inds = collect(indices(elem))
  v = CircularVector(inds)
  n = length(v)
  for i in 1:n
    half4pair[(v[i], v[i + 1])] = HalfEdge(v[i], eleminds[1])
  end

  # insert all other elements
  for e in 2:length(adjelems)
    elem = adjelems[e]
    inds = collect(indices(elem))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      # if pair of vertices is already in the
      # dictionary this means that the current
      # polygon has inconsistent orientation
      if haskey(half4pair, (v[i], v[i + 1]))
        # delete inserted pairs so far
        CCW[e] = false
        for j in 1:(i - 1)
          delete!(half4pair, (v[j], v[j + 1]))
        end
        break
      else
        # insert pair in consistent orientation
        half4pair[(v[i], v[i + 1])] = HalfEdge(v[i], eleminds[e])
      end
    end

    if !CCW[e]
      # reinsert pairs in CCW orientation
      for i in 1:n
        half4pair[(v[i + 1], v[i])] = HalfEdge(v[i + 1], eleminds[e])
      end
    end
  end

  # add missing pointers
  for (e, elem) in Iterators.enumerate(adjelems)
    inds = CCW[e] ? indices(elem) : reverse(indices(elem))
    v = CircularVector(collect(inds))
    n = length(v)
    for i in 1:n
      # update pointers prev and next
      he = half4pair[(v[i], v[i + 1])]
      he.prev = half4pair[(v[i - 1], v[i])]
      he.next = half4pair[(v[i + 1], v[i + 2])]

      # if not a border element, update half
      if haskey(half4pair, (v[i + 1], v[i]))
        he.half = half4pair[(v[i + 1], v[i])]
      else # create half-edge for border
        be = HalfEdge(v[i + 1], nothing)
        be.half = he
        he.half = be
      end
    end
  end

  # save halfedges in a vector of pairs
  halves = Vector{Tuple{HalfEdge,HalfEdge}}()
  visited = Set{Tuple{Int,Int}}()
  for ((u, v), he) in half4pair
    uv = minmax(u, v)
    if uv ∉ visited
      push!(halves, (he, he.half))
      push!(visited, uv)
    end
  end

  HalfEdgeTopology(halves, length(elems))
end

function adjsort(elems::AbstractVector{<:Connectivity})
  # initialize list of adjacent elements
  # with first element from original list
  list = indices.(elems)
  adjs = Tuple[popfirst!(list)]

  # the loop will terminate if the mesh
  # is manifold, and that is always true
  # with half-edge topology
  while !isempty(list)
    # lookup all elements that share at least
    # one vertex with the last adjacent element
    found = false
    vinds = last(adjs)
    for i in vinds
      einds = findall(e -> i ∈ e, list)
      if !isempty(einds)
        # lookup all elements that share at
        # least two vertices (i.e. edge)
        for j in sort(einds, rev=true)
          if length(vinds ∩ list[j]) > 1
            found = true
            push!(adjs, popat!(list, j))
          end
        end
      end
    end

    if !found && !isempty(list)
      # we are done with this connected component
      # pop a new element from the original list
      push!(adjs, popfirst!(list))
    end
  end

  connect.(adjs)
end

paramdim(::HalfEdgeTopology) = 2

"""
    half4elem(t, e)

Return a half-edge of the half-edge topology `t` on the `e`-th elem.
"""
half4elem(t::HalfEdgeTopology, e::Integer) = t.halfedges[t.half4elem[e]]

"""
    half4vert(t, v)

Return the half-edge of the half-edge topology `t` for which the
head is the `v`-th index.
"""
half4vert(t::HalfEdgeTopology, v::Integer) = t.halfedges[t.half4vert[v]]

"""
    half4edge(t, e)

Return the half-edge of the half-edge topology `t` for the edge `e`.

### Notes

Always return the half-edge to the "left".
"""
half4edge(t::HalfEdgeTopology, e::Integer) = t.halfedges[2e - 1]

"""
    half4pair(t, uv)

Return the half-edge of the half-edge topology `t` for the pair of
vertices `uv`.

### Notes

Always return the half-edge to the "left".
"""
half4pair(t::HalfEdgeTopology, uv::Tuple{Int,Int}) = half4edge(t, edge4pair(t, uv))

"""
    edge4pair(t, uv)

Return the edge of the half-edge topology `t` for the pair of vertices `uv`.
"""
edge4pair(t, uv) = t.edge4pair[minmax(uv...)]

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::HalfEdgeTopology) = length(t.half4vert)

function element(t::HalfEdgeTopology, ind)
  v = loop(half4elem(t, ind))
  connect(Tuple(v))
end

nelements(t::HalfEdgeTopology) = length(t.half4elem)

function facet(t::HalfEdgeTopology, ind)
  e = half4edge(t, ind)
  connect((e.head, e.half.head))
end

nfacets(t::HalfEdgeTopology) = length(t.halfedges) ÷ 2

# ------------
# CONVERSIONS
# ------------

Base.convert(::Type{HalfEdgeTopology}, t::Topology) = HalfEdgeTopology(collect(elements(t)))
