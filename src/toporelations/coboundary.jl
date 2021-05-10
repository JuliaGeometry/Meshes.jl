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
  𝒜 = Adjacency{0}(𝒞.structure)
  u, vs = vert, 𝒜(vert)
  elems = Int[]
  for v in vs
    e = half4pair((u, v), 𝒞.structure)
    h = e.half
    if e.elem ∉ elems
      push!(elems, e.elem)
    end
    if !isnothing(h.elem) && h.elem ∉ elems
      push!(elems, h.elem)
    end
  end
  elems
end

function (𝒞::Coboundary{1,2,S})(edge::Integer) where {S<:HalfEdgeStructure}
  e = half4edge(edge, 𝒞.structure)
  isnothing(e.half.elem) ? [e.elem] : [e.elem, e.half.elem]
end
