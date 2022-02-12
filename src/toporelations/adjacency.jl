# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Adjacency{P}(topology)

The adjacency relation of rank `P` for a given `topology`.
"""
struct Adjacency{P,D,T<:Topology} <: TopologicalRelation
  topology::T
end

function Adjacency{P}(topology) where {P}
  D = paramdim(topology)
  T = typeof(topology)

  @assert D ≥ P "invalid adjacency relation"

  Adjacency{P,D,T}(topology)
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (𝒜::Adjacency{0,2,T})(vert::Integer) where {T<:HalfEdgeTopology}
  e = half4vert(𝒜.topology, vert)

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
