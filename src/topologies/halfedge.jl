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
    HalfEdgeTopology(halfedges; nelems=nothing)

A data structure for orientable 2-manifolds based on
half-edges constructed from a vector of connectivity
`elements` or from a vector of pairs of `halfedges`.

The option `sort` can be used to sort the elements in
adjacent-first order in case of inconsistent orientation
(i.e. mix of clockwise and counter-clockwise).

The option `nelems` can be used to specify an approximate
number of `elements` as a size hint.

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

function HalfEdgeTopology(halves::AbstractVector{Tuple{HalfEdge,HalfEdge}}; nelems=nothing)
  # pre-allocate memory and provide size hints
  halfedges = Vector{HalfEdge}(undef, 2 * length(halves))
  edge4pair = Dict{Tuple{Int,Int},Int}()
  half4elem = Dict{Int,Int}()
  half4vert = Dict{Int,Int}()
  sizehint!(edge4pair, length(halves))
  if !isnothing(nelems)
    sizehint!(half4elem, nelems)
  end

  # flatten pairs of half-edges into a vector
  for (i, (h₁, h₂)) in enumerate(halves)
    # make sure that first half-edge is in the interior
    (h₁, h₂) = isnothing(h₁.elem) ? (h₂, h₁) : (h₁, h₂)

    # store half-edge with a linear index
    j = 2i - 1
    halfedges[j] = h₁
    halfedges[j + 1] = h₂

    # map element and vertex to a half-edge
    if !isnothing(h₁.elem)
      get!(half4elem, h₁.elem, j)
      get!(half4vert, h₁.head, j)
    end
    if !isnothing(h₂.elem)
      get!(half4elem, h₂.elem, j + 1)
      get!(half4vert, h₂.head, j + 1)
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
  elemsort = sort ? adjsort(elems) : elems
  eleminds = sort ? indexin(elemsort, elems) : 1:length(elems)
  adjelems::Vector{Vector{Int}} = map(collect ∘ indices, elemsort)

  # start assuming that all elements are
  # oriented consistently (e.g. CCW)
  isreversed = falses(length(adjelems))

  # initialize half-edges using first element
  half4pair = Dict{Tuple{Int,Int},HalfEdge}()
  elem = first(eleminds)
  inds = first(adjelems)
  n = length(inds)
  for i in eachindex(inds)
    u = inds[i]
    v = inds[mod1(i + 1, n)]
    # insert halfedges u -> v and v -> u
    he = get!(() -> HalfEdge(u, elem), half4pair, (u, v))
    half = get!(() -> HalfEdge(v, nothing), half4pair, (v, u))
    he.half = half
    half.half = he
  end

  # update half-edges using all other elements
  added = false
  disconnected = false
  remaining = collect(2:length(adjelems))
  while !isempty(remaining)
    iter = 1
    while iter ≤ length(remaining)
      other = remaining[iter]
      elem = eleminds[other]
      inds = adjelems[other]
      n = length(inds)
      if anyhalf(inds, half4pair) || disconnected
        # at least one half-edge has been found, so we can assess
        # the orientation w.r.t. previously added elements/edges
        added = true
        disconnected = false
        deleteat!(remaining, iter)
      else
        iter += 1
        continue
      end

      isreversed[other] = anyhalfclaimed(inds, half4pair)

      # insert half-edges in consistent orientation
      if isreversed[other]
        for i in eachindex(inds)
          u = inds[i]
          v = inds[mod1(i + 1, n)]
          he = get!(() -> HalfEdge(v, elem), half4pair, (v, u))
          if isnothing(he.elem)
            he.elem = elem
          else
            assertion(he.elem === elem, lazy"inconsistent duplicate edge $he for $(elem) and $(he.elem)")
          end
          half = get!(() -> HalfEdge(u, nothing), half4pair, (u, v))
          he.half = half
          half.half = he
        end
      else
        for i in eachindex(inds)
          u = inds[i]
          v = inds[mod1(i + 1, n)]
          he = get!(() -> HalfEdge(u, elem), half4pair, (u, v))
          he.elem = elem # this may be a pre-allocated half-edge with a nothing `elem`
          half = get!(() -> HalfEdge(v, nothing), half4pair, (v, u))
          he.half = half
          half.half = he
        end
      end
    end

    if added
      added = false
    elseif !isempty(remaining)
      disconnected = true
      added = false
    end
  end

  # add missing pointers and save halfedges in a vector of pairs
  halves = Vector{Tuple{HalfEdge,HalfEdge}}()
  visited = Set{Tuple{Int,Int}}()
  for (e, inds) in enumerate(adjelems)
    inds = isreversed[e] ? circshift!(reverse!(inds), 1) : inds
    n = length(inds)
    for i in eachindex(inds)
      u = inds[i]
      v = inds[mod1(i + 1, n)]
      w = inds[mod1(i + 2, n)]

      # update pointers prev and next
      he = half4pair[(u, v)]
      he.next = half4pair[(v, w)]
      he.next.prev = he

      uv = minmax(u, v)
      if uv ∉ visited
        push!(halves, (he, he.half))
        push!(visited, uv)
      end
    end
  end

  HalfEdgeTopology(halves; nelems=length(elems))
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

# -----------------
# HELPER FUNCTIONS
# -----------------

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

function anyhalf(inds, half4pair)
  n = length(inds)
  for i in eachindex(inds)
    uv = (inds[i], inds[mod1(i + 1, n)])
    if haskey(half4pair, uv)
      return true
    end
  end
  return false
end

function anyhalfclaimed(inds, half4pair)
  n = length(inds)
  for i in eachindex(inds)
    uv = (inds[i], inds[mod1(i + 1, n)])
    if !isnothing(get(() -> HalfEdge(0, nothing), half4pair, uv).elem)
      return true
    end
  end
  return false
end
