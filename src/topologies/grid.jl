# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GridTopology(dims, [open])

A data structure for grid topologies with `dims` elements.
Optionally, specify which dimensions are `open` and which
are closed, i.e. do not wrap around.

## Examples

```julia
julia> GridTopology((10,20)) # 10x20 elements in a grid
julia> GridTopology((10,20), (true,false)) # cylinder topology
julia> GridTopology((10,20), (false,false)) # sphere topology
```
"""
struct GridTopology{D} <: Topology
  dims::Dims{D}
  open::NTuple{D,Bool}
end

GridTopology(dims::Dims{D}) where {D} = GridTopology{D}(dims, ntuple(i->true, D))

GridTopology(dims::Vararg{Int,D}) where {D} = GridTopology(dims)

paramdim(::GridTopology{D}) where {D} = D

Base.size(t::GridTopology) = t.dims

"""
    isclosed(t)

Tells whether or not the grid topology `t` is closed
along each dimension.
"""
isclosed(t::GridTopology) = .!t.open

"""
    elem2cart(t, e)

Return the Cartesian indices of the `e`-th element of
the grid topology `t`.
"""
elem2cart(t::GridTopology, e) = CartesianIndices(t.dims)[e].I

"""
    cart2elem(t, i, j, k, ...)

Return the linear index of the element of the
grid topology `t` with Cartesian indices
`i`, `j`, `k`, ...
"""
cart2elem(t::GridTopology, ijk...) = LinearIndices(t.dims)[ijk...]

"""
    corner2cart(t, v)

Return the Cartesian indices of the element of the
grid topology `t` with top left vertex `v`.
"""
corner2cart(t, v) = CartesianIndices(t.dims .+ t.open)[v].I

"""
    cart2corner(t, i, j, k, ...)

Return the linear index of the top left vertex of the
element of the grid topology `t` with Cartesian indices
`i`, `j`, `k`, ...
"""
cart2corner(t::GridTopology, ijk...) = LinearIndices(t.dims .+ t.open)[ijk...]

"""
    elem2corner(t, e)

Return the linear index of the top left vertex of the
`e`-th element of the grid topology `t`.
"""
elem2corner(t::GridTopology, e) = cart2corner(t, elem2cart(t, e)...)

"""
    corner2elem(t, v)

Return the linear index of the element of the grid
topology `t` with top left vertex `v`
"""
corner2elem(t::GridTopology, v) = cart2elem(t, corner2cart(t, v)...)

"""
    rank2type(t, rank)

Return the polytope type of given `rank` for the grid
topology `t`.
"""
function rank2type(::GridTopology{D}, rank::Integer) where {D}
  @assert rank ≤ D "invalid rank for grid topology"
  rank == 1 && return Segment
  rank == 2 && return Quadrangle
  rank == 3 && return Hexahedron
  throw(ErrorException("not implemented"))
end

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::GridTopology) = prod(t.dims .+ t.open)

function faces(t::GridTopology{D}, rank) where {D}
  if rank == D
    elements(t)
  elseif rank == D - 1
    facets(t)
  else
    throw(ErrorException("not implemented"))
  end
end

function element(t::GridTopology{D}, ind) where {D}
  ∂ = Boundary{D,0}(t)
  T = rank2type(t, D)
  connect(Tuple(∂(ind)), T)
end

nelements(t::GridTopology) = prod(t.dims)

function facet(t::GridTopology{D}, ind) where {D}
  ∂ = Boundary{D-1,0}(t)
  T = rank2type(t, D-1)
  connect(Tuple(∂(ind)), T)
end

nfacets(t::GridTopology{1}) =
  t.dims[1] + t.open[1]

nfacets(t::GridTopology{2}) =
  2prod(t.dims) +
  t.open[2]*t.dims[1] +
  t.open[1]*t.dims[2]

nfacets(t::GridTopology{3}) =
  3prod(t.dims) +
  t.open[3]*(t.dims[1]*t.dims[2]) +
  t.open[2]*(t.dims[1]*t.dims[3]) +
  t.open[1]*(t.dims[2]*t.dims[3])