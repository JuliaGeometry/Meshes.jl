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
    HalfEdgeStructure(halfedges, edgeonelem, edgeonvert)

A data structure for orientable 2-manifolds based
on half-edges.

Two types of half-edges exist (Kettner 1999). This
implementation is the most common type that splits
the incident elements.

A vector of `halfedges` together with a dictionary of
`edgeonelem` and a dictionary of `edgeonvert` can be
used to retrieve topolological relations in optimal
time. In this case, `edgeonvert[i]` returns the index
of the half-edge in `halfedges` with head equal to `i`.
Similarly, `edgeonelem[i]` returns the index of a
half-edge in `halfedges` that is in the element `i`.

See also [`TopologicalStructure`](@ref).

## References

* Kettner, L. (1999). [Using generic programming for
  designing a data structure for polyhedral surfaces]
  (https://www.sciencedirect.com/science/article/pii/S0925772199000073)
"""
struct HalfEdgeStructure <: TopologicalStructure
  halfedges::Vector{HalfEdge}
  edgeonelem::Dict{Int,Int}
  edgeonvert::Dict{Int,Int}
end

function HalfEdgeStructure(elems::AbstractVector{<:Connectivity})
  @assert all(e -> paramdim(e) == 2, elems) "invalid element for half-edge structure"

  # number of vertices is the maximum index in connectivities
  nvertices = maximum(i for e in elems for i in indices(e))

  # initialization step
  edge4pair = Dict{Tuple{Int,Int},HalfEdge}()
  for (e, elem) in Iterators.enumerate(elems)
    inds = collect(indices(elem))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      edge4pair[(v[i], v[i+1])] = HalfEdge(v[i], e)
    end
  end

  # add missing pointers
  for elem in elems
    inds = collect(indices(elem))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      # update pointers prev and next
      he = edge4pair[(v[i], v[i+1])]
      he.prev = edge4pair[(v[i-1],   v[i])]
      he.next = edge4pair[(v[i+1], v[i+2])]

      # if not a border element, update half
      if haskey(edge4pair, (v[i+1], v[i]))
        he.half = edge4pair[(v[i+1], v[i])]
      else # create half-edge for border
        be = HalfEdge(v[i+1], nothing)
        be.half = he
        he.half = be
      end
    end
  end

  # store all halfedges in a vector
  halfedges = HalfEdge[]
  processed = Set{Tuple{Int,Int}}()
  for ((u, v), he) in edge4pair
    if (u, v) ∉ processed
      append!(halfedges, [he, he.half])
      push!(processed, (u, v))
      push!(processed, (v, u))
    end
  end

  # reverse mappings
  edgeonelem = Dict{Int,Int}()
  edgeonvert = Dict{Int,Int}()
  for (e, he) in enumerate(halfedges)
    if !isnothing(he.elem) # interior half-edge
      if !haskey(edgeonelem, he.elem)
        edgeonelem[he.elem] = e
      end
      if !haskey(edgeonvert, he.head)
        edgeonvert[he.head] = e
      end
    end
  end

  HalfEdgeStructure(halfedges, edgeonelem, edgeonvert)
end

"""
    edgeonelem(e, s)

Return a half-edge of the half-edge structure `s` on the `e`-th elem.
"""
edgeonelem(e::Integer, s::HalfEdgeStructure) = s.halfedges[s.edgeonelem[e]]

"""
    edgeonvert(v, s)

Return the half-edge of the half-edge structure `s` for which the
head is the `v`-th index.
"""
edgeonvert(v::Integer, s::HalfEdgeStructure) = s.halfedges[s.edgeonvert[v]]

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
  eᵤ = edgeonvert(u, s)

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
    [ngon4edge(e.half)]
  elseif isnothing(e.half.elem)
    [ngon4edge(e)]
  else
    [ngon4edge(e), ngon4edge(e.half)]
  end
end

function adjacency(c::Connectivity{<:Polygon}, s::HalfEdgeStructure)
end

function adjacency(c::Connectivity{<:Segment}, s::HalfEdgeStructure)
end

function adjacency(v::Integer, s::HalfEdgeStructure)
  e = edgeonvert(v, s)
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

# helper function to construct
# n-gon from given half-edge
function ngon4edge(e::HalfEdge)
  n = e.next
  v = [e.head]
  while n != e
    push!(v, n.head)
    n = n.next
  end
  connect(Tuple(v), Ngon{length(v)})
end

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

element(s::HalfEdgeStructure, ind) = ngon4edge(edgeonelem(ind, s))

nelements(s::HalfEdgeStructure) = length(s.edgeonelem)

function facet(s::HalfEdgeStructure, ind)
  e = s.halfedges[2ind-1]
  connect((e.head, e.half.head), Segment)
end

nfacets(s::HalfEdgeStructure) = length(s.halfedges) ÷ 2

# ---------------------------------
# CONVERSION FROM OTHER STRUCTURES
# ---------------------------------

Base.convert(::Type{HalfEdgeStructure}, s::TopologicalStructure) =
  HalfEdgeStructure(collect(elements(s)))
