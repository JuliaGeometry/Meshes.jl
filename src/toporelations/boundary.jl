# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Boundary{P,Q}(topology)

The boundary relation from rank `P` to smaller rank `Q` for
a given `topology`.
"""
struct Boundary{P,Q,D,T<:Topology} <: TopologicalRelation
  topology::T
end

function Boundary{P,Q}(topology) where {P,Q}
  D = paramdim(topology)
  T = typeof(topology)

  @assert D ≥ P > Q "invalid boundary relation"

  Boundary{P,Q,D,T}(topology)
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (∂::Boundary{2,1,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  t = ∂.topology
  l = loop(half4elem(elem, t))
  v = CircularVector(l)
  [edge4pair((v[i], v[i+1]), t) for i in 1:length(v)]
end

function (∂::Boundary{2,0,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  loop(half4elem(elem, ∂.topology))
end

function (∂::Boundary{1,0,2,T})(edge::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(edge, ∂.topology)
  [e.head, e.half.head]
end
