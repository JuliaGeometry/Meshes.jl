# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Coboundary{P,Q,S}

The co-boundary relation from rank `P` to greater rank `Q` for
topological structure of type `S`.
"""
struct Coboundary{P,Q,S<:TopologicalStructure} <: TopologicalRelation
  structure::S
end

Coboundary{P,Q}(structure::S) where {P,Q,S} = Coboundary{P,Q,S}(structure)

# --------------------
# HALF-EDGE STRUCTURE
# --------------------

function (𝒞::Coboundary{0,1,S})(vert::Integer) where {S<:HalfEdgeStructure}
  s = 𝒞.structure
  𝒜 = Adjacency{0}(s)
  [edge4pair((vert, other), s) for other in 𝒜(vert)]
end

function (𝒞::Coboundary{0,2,S})(vert::Integer) where {S<:HalfEdgeStructure}
  e = half4vert(vert, 𝒞.structure)
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
  o = p.half
  elems = [e.elem]
  while !isnothing(o.elem) && o != e
    push!(elems, o.elem)
    p = o.prev
    o = p.half
  end

  elems
end

function (𝒞::Coboundary{1,2,S})(edge::Integer) where {S<:HalfEdgeStructure}
  e = half4edge(edge, 𝒞.structure)
  isnothing(e.half.elem) ? [e.elem] : [e.elem, e.half.elem]
end
