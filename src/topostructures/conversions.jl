# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Base.convert(::Type{HalfEdgeStructure}, s::ExplicitStructure)
  # half-edge structure only works with orientable 2-manifolds
  elems = (c for c in connectivities(s) if paramdim(c) == 2)
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
  inneredges  = HalfEdge[]
  borderedges = HalfEdge[]
  edgeonelem  = Int[]
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
        push!(borderedges, be)
      end

      # save inner half-edges in order
      push!(inneredges, he)
    end

    # save first halfedge for this elem
    push!(edgeonelem, length(inneredges) - n + 1)
  end

  # all half-edges (inner and border)
  halfedges = [inneredges; borderedges]

  # map vertices to inner half-edges
  edgeonvertex = Vector{Int}(undef, nvertices)
  for i in eachindex(edgeonvertex)
    e = findfirst(he -> he.head == i, inneredges)
    edgeonvertex[i] = e
  end

  HalfEdgeStructure(halfedges, edgeonelem, edgeonvertex)
end

function Base.convert(::Type{ExplicitStructure}, s::HalfEdgeStructure)
  connec = map(1:nelements(s)) do i
    # select a half-edge on the elem
    e = edgeonelem(s, i)

    # retrieve vertices of the element
    n = e.next
    v = [e.head]
    while n != e
      push!(v, n.head)
      n = n.next
    end

    # connect vertices into a n-gon
    connect(Tuple(v), Ngon{length(v)})
  end

  ExplicitStructure(connec)
end
