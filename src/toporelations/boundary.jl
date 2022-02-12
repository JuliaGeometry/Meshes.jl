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

# --------------
# GRID TOPOLOGY
# --------------

function (∂::Boundary{3,0,3,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  i, j, k = elem2cart(t, ind)
  i1 = cart2corner(t, i  , j  , k  )
  i2 = cart2corner(t, i+1, j  , k  )
  i3 = cart2corner(t, i+1, j+1, k  )
  i4 = cart2corner(t, i  , j+1, k  )
  i5 = cart2corner(t, i  , j  , k+1)
  i6 = cart2corner(t, i+1, j  , k+1)
  i7 = cart2corner(t, i+1, j+1, k+1)
  i8 = cart2corner(t, i  , j+1, k+1)
  [i1, i2, i3, i4, i5, i6, i7, i8]
end

function (∂::Boundary{2,0,2,T})(ind::Integer) where {T<:GridTopology}
  t = ∂.topology
  i, j = elem2cart(t, ind)
  i1 = cart2corner(t, i  , j  )
  i2 = cart2corner(t, i+1, j  )
  i3 = cart2corner(t, i+1, j+1)
  i4 = cart2corner(t, i  , j+1)
  [i1, i2, i3, i4]
end

function (∂::Boundary{1,0,1,T})(ind::Integer) where {T<:GridTopology}
  i1 = ind
  i2 = ind+1
  [i1, i2]
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

function (∂::Boundary{2,1,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  t = ∂.topology
  l = loop(half4elem(t, elem))
  v = CircularVector(l)
  [edge4pair(t, (v[i], v[i+1])) for i in 1:length(v)]
end

function (∂::Boundary{2,0,2,T})(elem::Integer) where {T<:HalfEdgeTopology}
  loop(half4elem(∂.topology, elem))
end

function (∂::Boundary{1,0,2,T})(edge::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(∂.topology, edge)
  [e.head, e.half.head]
end
