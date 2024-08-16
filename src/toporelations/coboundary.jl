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

  inds = Int[]
  for offset in CartesianIndices(ntuple(i -> -1:0, D))
    # apply offset to center index
    sind = cind .+ Tuple(offset)

    # wrap indices in case of periodic dimension
    wind = ntuple(D) do i
      cycl[i] ? mod1(sind[i], dims[i]) : sind[i]
    end

    # discard invalid indices
    valid(i) = 1 ≤ wind[i] ≤ dims[i]
    all(valid, 1:D) && push!(inds, cart2elem(topo, wind...))
  end

  ntuple(i -> inds[i], length(inds))
end

# -------------------
# HALF-EDGE TOPOLOGY
# -------------------

# segments sharing a vertex in 2D mesh
function (𝒞::Coboundary{0,1,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  t = 𝒞.topology
  𝒜 = Adjacency{0}(t)
  o = 𝒜(ind)
  ntuple(length(o)) do i
    edge4pair(t, (ind, o[i]))
  end
end

# elements sharing a vertex in 2D mesh
function (𝒞::Coboundary{0,2,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  e = half4vert(𝒞.topology, ind)

  # initialize result
  inds = [e.elem]

  # search in CCW orientation
  p = e.prev
  h = p.half
  while !isnothing(h.elem) && h != e
    push!(inds, h.elem)
    p = h.prev
    h = p.half
  end

  # if border edge is hit
  if isnothing(h.elem)
    # search in CW orientation
    h = e.half
    while !isnothing(h.elem)
      pushfirst!(inds, h.elem)
      n = h.next
      h = n.half
    end
  end

  ntuple(i -> inds[i], length(inds))
end

# elements sharing a segment in 2D mesh
function (𝒞::Coboundary{1,2,2,T})(ind::Integer) where {T<:HalfEdgeTopology}
  e = half4edge(𝒞.topology, ind)
  isnothing(e.half.elem) ? (e.elem,) : (e.elem, e.half.elem)
end
