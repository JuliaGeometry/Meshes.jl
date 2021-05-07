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
  𝒜 = Adjacency{0}(𝒞.structure)
  [(vert, other) for other in 𝒜(vert)]
end

function (𝒞::Coboundary{0,2,S})(vert::Integer) where {S<:HalfEdgeStructure}
end

function (𝒞::Coboundary{1,2,S})(edge::Integer) where {S<:HalfEdgeStructure}
  e = half4edge(edge, 𝒞.structure)
  elems = if isnothing(e.elem)
    [loop(e.half)]
  elseif isnothing(e.half.elem)
    [loop(e)]
  else
    [loop(e), loop(e.half)]
  end
  Tuple.(elems)
end
