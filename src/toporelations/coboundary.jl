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

  assertion(P < Q ≤ D, "invalid coboundary relation")

  Coboundary{P,Q,D,T}(topology)
end

# --------------
# GRID TOPOLOGY
# --------------

# elements sharing vertex in grid
function (𝒞::Coboundary{0,D,D,T})(ind::Integer) where {D,T<:GridTopology}
  topo = 𝒞.topology
  dims = size(topo)
  cycl = isperiodic(topo)
  cind = corner2cart(topo, ind)

  # offsets along each dimension
  offsets = CartesianIndices(ntuple(i -> -1:0, D))

  ninds = NTuple{D,Int}[]
  for offset in offsets
    # apply offset to center index
    sind = cind .+ Tuple(offset)

    # wrap indices in case of periodic dimension
    wrap(i) = mod1(sind[i], dims[i])
    wind = ntuple(i -> cycl[i] ? wrap(i) : sind[i], D)

    # discard invalid indices
    valid(i) = 1 ≤ wind[i] ≤ dims[i]
    all(valid, 1:D) && push!(ninds, wind)
  end

  # return linear index of element
  [cart2elem(topo, ind...) for ind in ninds]
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

# segments sharing a vertex in 2D mesh
function (𝒞::Coboundary{0,1,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  t = 𝒞.topology
  𝒜 = Adjacency{0}(t)
  [edge4pair(t, (ind, other)) for other in 𝒜(ind)]
end

# elements sharing a vertex in 2D mesh
function (𝒞::Coboundary{0,2,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  e = half4vert(𝒞.topology, ind)

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

# elements sharing a segment in 2D mesh
function (𝒞::Coboundary{1,2,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(𝒞.topology, ind)
  isnothing(e.half.elem) ? [e.elem] : [e.elem, e.half.elem]
end
