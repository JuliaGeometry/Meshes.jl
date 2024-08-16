# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GridTopology(dims, [periodic])

A data structure for grid topologies with `dims` elements.
Optionally, specify which dimensions are `periodic`. Default
to aperiodic dimensions.

## Examples

```julia
julia> GridTopology((10,20)) # 10x20 elements in a grid
julia> GridTopology((10,20), (true,false)) # cylinder topology
```
"""
struct GridTopology{D} <: Topology
  dims::Dims{D}
  open::NTuple{D,Bool}

  function GridTopology{D}(dims, periodic) where {D}
    new(dims, .!periodic)
  end
end

GridTopology(dims, periodic) = GridTopology{length(dims)}(dims, periodic)

GridTopology(dims::Dims{D}) where {D} = GridTopology(dims, ntuple(i -> false, D))

GridTopology(dims::Vararg{Int,D}) where {D} = GridTopology(dims)

paramdim(::GridTopology{D}) where {D} = D

Base.size(t::GridTopology) = t.dims

"""
    isperiodic(topology)

Tells whether or not the `topology` is periodic
along each parametric dimension.
"""
isperiodic(t::GridTopology) = .!t.open

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
corner2cart(t::GridTopology, v) = CartesianIndices(t.dims .+ t.open)[v].I

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
    elementtype(t)

Return the element type of the grid topology `t`.
"""
elementtype(::GridTopology{1}) = Segment
elementtype(::GridTopology{2}) = Quadrangle
elementtype(::GridTopology{3}) = Hexahedron

"""
    facettype(t)

Return the facet type of the grid topology `t`.
"""
facettype(::GridTopology{1}) = Point
facettype(::GridTopology{2}) = Segment
facettype(::GridTopology{3}) = Quadrangle

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::GridTopology) = prod(t.dims .+ t.open)

function element(t::GridTopology{D}, ind) where {D}
  ∂ = Boundary{D,0}(t)
  T = elementtype(t)
  connect(∂(ind), T)
end

nelements(t::GridTopology) = prod(t.dims)

function facet(t::GridTopology{D}, ind) where {D}
  ∂ = Boundary{D - 1,0}(t)
  T = facettype(t)
  connect(∂(ind), T)
end

nfacets(t::GridTopology{1}) = t.dims[1] + t.open[1]

nfacets(t::GridTopology{2}) = 2prod(t.dims) + t.open[2] * t.dims[1] + t.open[1] * t.dims[2]

nfacets(t::GridTopology{3}) =
  3prod(t.dims) +
  t.open[3] * (t.dims[1] * t.dims[2]) +
  t.open[2] * (t.dims[1] * t.dims[3]) +
  t.open[1] * (t.dims[2] * t.dims[3])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, t::GridTopology)
  dims = join(t.dims, "×")
  strs = replace(t.open, true => "aperiodic", false => "periodic")
  peri = join(strs, ", ")
  print(io, "$dims GridTopology($peri)")
end
