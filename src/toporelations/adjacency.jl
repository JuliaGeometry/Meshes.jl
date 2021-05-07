# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Adjacency{P,S}

The adjacency relation of rank `P` for topological structure of type `S`.
"""
struct Adjacency{P,S<:TopologicalStructure} <: TopologicalRelation
  structure::S
end

Adjacency{P}(structure::S) where {P,S} = Adjacency{P,S}(structure)

# --------------------
# HALF-EDGE STRUCTURE
# --------------------

function (𝒜::Adjacency{0,S})(vert::Integer) where {S<:HalfEdgeStructure}
  e = half4vert(vert, 𝒜.structure)
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
  v = [n.head]
  while !isnothing(o.elem) && o != e
    p = o.prev
    n = o.next
    o = p.half
    push!(v, n.head)
  end
  # if border edge is hit, add last arm manually
  isnothing(o.elem) && push!(v, o.half.head)

  v
end
