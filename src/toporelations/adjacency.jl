# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Adjacency{P,T}

The adjacency relation of rank `P` for topology of type `T`.
"""
struct Adjacency{P,T<:Topology} <: TopologicalRelation
  topology::T
end

Adjacency{P}(topology::T) where {P,T} = Adjacency{P,T}(topology)

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (𝒜::Adjacency{0,T})(vert::Integer) where {T<:HalfEdgeTopology}
  e = half4vert(vert, 𝒜.topology)

  # initialize result
  vertices = [e.half.head]

  # search in CCW orientation
  p = e.prev
  h = p.half
  while !isnothing(h.elem) && h != e
    push!(vertices, p.head)
    p = h.prev
    h = p.half
  end

  # if border edge is hit
  if isnothing(h.elem)
    # add last arm manually
    push!(vertices, p.head)

    # search in CW orientation
    h = e.half
    while !isnothing(h.elem)
      n = h.next
      h = n.half
      pushfirst!(vertices, h.head)
    end
  end

  vertices
end
