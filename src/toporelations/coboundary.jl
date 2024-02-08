# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Coboundary{P,Q}(topology)

The co-boundary relation from rank `P` to greater rank `Q` for
a given `topology`.
"""
struct Coboundary{P,Q,D,T<:Topology} <: TopologicalRelation
  topology::T
end

function Coboundary{P,Q}(topology) where {P,Q}
  D = paramdim(topology)
  T = typeof(topology)

  @assert P < Q ≤ D "invalid coboundary relation"

  Coboundary{P,Q,D,T}(topology)
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (𝒞::Coboundary{0,1,2,T})(vert::Integer) where {T<:HalfEdgeTopology}
  t = 𝒞.topology
  𝒜 = Adjacency{0}(t)
  [edge4pair(t, (vert, other)) for other in 𝒜(vert)]
end

function (𝒞::Coboundary{0,2,2,T})(vert::Integer) where {T<:HalfEdgeTopology}
  e = half4vert(𝒞.topology, vert)

  # initialize result
  elements = [e.elem]

  # search in CCW orientation
  p = e.prev
  h = p.half
  while !isnothing(h.elem) && h != e
    push!(elements, h.elem)
    p = h.prev
    h = p.half
  end

  # if border edge is hit
  if isnothing(h.elem)
    # search in CW orientation
    h = e.half
    while !isnothing(h.elem)
      pushfirst!(elements, h.elem)
      n = h.next
      h = n.half
    end
  end

  elements
end

function (𝒞::Coboundary{1,2,2,T})(edge::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(𝒞.topology, edge)
  isnothing(e.half.elem) ? [e.elem] : [e.elem, e.half.elem]
end

# -------------------
# IndexedAdjacenciesTopology
# -------------------

function (𝒞::Coboundary{0,K,K,T})(vert_idx::Int) where {K, T<:IndexedAdjacenciesTopology}
    𝒜 = Adjacency{K}(𝒞.topology)
    results = Int[]
    indices_to_explore = Int[]  # we treat this like a queue
    push!(indices_to_explore, 𝒞.topology.R_star_relations[vert_idx])
    while !isempty(indices_to_explore)
        idx = popfirst!(indices_to_explore)
        if idx == -1 || idx ∈ results
            continue
        elseif vert_idx ∈ indices(𝒞.topology.simplicies[idx])
            push!(results, idx)
            append!(indices_to_explore, 𝒜(idx))
        end
    end
    results
end
