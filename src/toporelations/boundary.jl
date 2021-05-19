# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Boundary{P,Q,T}

The boundary relation from rank `P` to smaller rank `Q` for
topology of type `T`.
"""
struct Boundary{P,Q,T<:Topology} <: TopologicalRelation
  topology::T
end

Boundary{P,Q}(topology::T) where {P,Q,T} = Boundary{P,Q,T}(topology)

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (∂::Boundary{2,1,T})(elem::Integer) where {T<:HalfEdgeTopology}
  t = ∂.topology
  l = loop(half4elem(elem, t))
  v = CircularVector(l)
  [edge4pair((v[i], v[i+1]), t) for i in 1:length(v)]
end

function (∂::Boundary{2,0,T})(elem::Integer) where {T<:HalfEdgeTopology}
  loop(half4elem(elem, ∂.topology))
end

function (∂::Boundary{1,0,T})(edge::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(edge, ∂.topology)
  [e.head, e.half.head]
end
