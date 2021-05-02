# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Base.convert(::Type{HalfEdgeStructure}, s::ElementListStructure)
  # half-edge structure only works with orientable 2-manifolds
  elems = elements(s)
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

function Base.convert(::Type{ElementListStructure}, s::HalfEdgeStructure)
  ElementListStructure(collect(elements(s)))
end
