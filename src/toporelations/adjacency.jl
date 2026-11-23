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

  assertion(D ≥ P, "invalid adjacency relation")

  Adjacency{P,D,T}(topology)
end

# --------------
# GRID TOPOLOGY
# --------------

# adjacent vertices in a D-dimensional grid topology
function (𝒜::Adjacency{0,D,T})(ind::Int) where {D,T<:GridTopology}
  # retrieve topology info
  topo = 𝒜.topology
  dims = size(topo)
  cycl = isperiodic(topo)

  # construct topology for vertices
  vtopo = GridTopology(dims .+ 1, cycl)
  𝒜vert = Adjacency{D}(vtopo)

  𝒜vert(ind)
end

# adjacent elements in a D-dimensional grid topology
function (𝒜::Adjacency{D,D,T})(ind::Int) where {D,T<:GridTopology}
  topo = 𝒜.topology
  dims = size(topo)
  cycl = isperiodic(topo)
  cind = elem2cart(topo, ind)

  inds = Int[]
  for d in 1:D, s in (-1, 1)
    # offset along each dimension
    offset = ntuple(i -> i == d ? s : 0, D)

    # apply offset to center index
    sind = cind .+ offset

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

# adjacent vertices in a 2D half-edge topology
function (𝒜::Adjacency{0,2,T})(vert::Int) where {T<:HalfEdgeTopology}
  e = half4vert(𝒜.topology, vert)

  # initialize result
  inds = [head(half(e))]

  # search in CCW orientation
  p = prev(e)
  h = half(p)
  while !isnothing(elem(h)) && h ≠ e
    push!(inds, head(p))
    p = prev(h)
    h = half(p)
  end

  # if border edge is hit
  if isnothing(elem(h))
    # add last arm manually
    push!(inds, head(p))

    # search in CW orientation
    h = half(e)
    while !isnothing(elem(h))
      n = next(h)
      h = half(n)
      pushfirst!(inds, head(h))
    end
  end

  ntuple(i -> inds[i], length(inds))
end

# adjacent elements in a 2D half-edge topology
function (𝒜::Adjacency{2,2,T})(ind::Int) where {T<:HalfEdgeTopology}
  inds = Int[]

  e = half4elem(𝒜.topology, ind)
  i = elem(half(e))
  isnothing(i) || push!(inds, i)

  n = next(e)
  while n ≠ e
    i = elem(half(n))
    isnothing(i) || push!(inds, i)
    n = next(n)
  end

  ntuple(i -> inds[i], length(inds))
end
