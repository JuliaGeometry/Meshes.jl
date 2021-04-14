# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Base.convert(::Type{HalfEdgeStructure}, s::ExplicitStructure)
  # half-edge structure only works with orientable 2-manifolds
  cells = (c for c in connectivities(s) if paramdim(c) == 2)
  nvertices = maximum(i for c in cells for i in indices(c))

  # initialization step
  edge4pair = Dict{Tuple{Int,Int},HalfEdge}()
  for (c, cell) in Iterators.enumerate(cells)
    inds = collect(indices(cell))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      edge4pair[(v[i], v[i+1])] = HalfEdge(v[i], c)
    end
  end

  # add missing pointers
  inneredges  = HalfEdge[]
  borderedges = HalfEdge[]
  edgeoncell  = Int[]
  for (c, cell) in Iterators.enumerate(cells)
    inds = collect(indices(cell))
    v = CircularVector(inds)
    n = length(v)
    for i in 1:n
      # update pointers prev and next
      he = edge4pair[(v[i], v[i+1])]
      he.prev = edge4pair[(v[i-1],   v[i])]
      he.next = edge4pair[(v[i+1], v[i+2])]

      # if not a border cell, update half
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

    # save first halfedge for this cell
    push!(edgeoncell, length(inneredges) - n + 1)
  end

  # all half-edges (inner and border)
  halfedges = [inneredges; borderedges]

  # map vertices to inner half-edges
  edgeonvertex = Vector{Int}(undef, nvertices)
  for i in eachindex(edgeonvertex)
    e = findfirst(he -> he.head == i, inneredges)
    edgeonvertex[i] = e
  end

  HalfEdgeStructure(halfedges, edgeoncell, edgeonvertex)
end

function Base.convert(::Type{ExplicitStructure}, s::HalfEdgeStructure)
  connec = map(1:ncells(s)) do c
    # select a half-edge on the cell
    e = edgeoncell(s, c)

    # retrieve vertices of the cell
    n = e.next
    v = [e.head]
    while n != e
      push!(v, n.head)
      n = n.next
    end

    # connect vertices into a polytope type
    nv = length(v)
    if nv == 3
      celltype = Triangle
    elseif nv == 4
      celltype = Quadrangle
    else
      throw(Error("Polytope type not implemented."))
    end

    connect(Tuple(v), celltype)
  end

  ExplicitStructure(connec)
end
