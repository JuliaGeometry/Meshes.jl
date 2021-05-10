# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Boundary{P,Q,S}

The boundary relation from rank `P` to smaller rank `Q` for
topological structure of type `S`.
"""
struct Boundary{P,Q,S<:TopologicalStructure} <: TopologicalRelation
  structure::S
end

Boundary{P,Q}(structure::S) where {P,Q,S} = Boundary{P,Q,S}(structure)

# --------------------
# HALF-EDGE STRUCTURE
# --------------------

function (∂::Boundary{2,1,S})(elem::Integer) where {S<:HalfEdgeStructure}
  s = ∂.structure
  l = loop(half4elem(elem, s))
  v = CircularVector(l)
  [edge4pair((v[i], v[i+1]), s) for i in 1:length(v)]
end

function (∂::Boundary{2,0,S})(elem::Integer) where {S<:HalfEdgeStructure}
  loop(half4elem(elem, ∂.structure))
end
