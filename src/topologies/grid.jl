# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GridTopology(dims)

A data structure for grid topologies with `dims` elements.
"""
struct GridTopology{D} <: Topology
  dims::Dims{D}
end

GridTopology(dims::Vararg{Int,D}) where {D} = GridTopology{D}(dims)

paramdim(::GridTopology{D}) where {D} = D

"""
    size(t)

Return the size of the grid topology `t`, i.e. the
number of elements along each dimension of the grid.
"""
Base.size(t::GridTopology) = t.dims

"""
    i, j, k, ... = elem2cart(t, e)

Return the Cartesian indices of the `e`-th element of
the grid topology `t`.
"""
elem2cart(t::GridTopology, e) = CartesianIndices(t.dims)[e].I

"""
    cart2corner(t, i, j, k, ...)

Return the linear index of the top left vertex of the
element of the topology `t` with Cartesian indices
`i`, `j`, `k`, ...
"""
cart2corner(t::GridTopology, ijk...) = LinearIndices(t.dims .+ 1)[ijk...]

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::GridTopology) = prod(t.dims .+ 1)

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
  if D == 1
    i1 = ind
    i2 = ind+1
    connect((i1, i2), Segment)
  elseif D == 2
    i, j = elem2cart(t, ind)
    i1 = cart2corner(t, i  , j  )
    i2 = cart2corner(t, i+1, j  )
    i3 = cart2corner(t, i+1, j+1)
    i4 = cart2corner(t, i  , j+1)
    connect((i1, i2, i3, i4), Quadrangle)
  elseif D == 3
    i, j, k = elem2cart(t, ind)
    i1 = cart2corner(t, i  , j  , k  )
    i2 = cart2corner(t, i+1, j  , k  )
    i3 = cart2corner(t, i+1, j+1, k  )
    i4 = cart2corner(t, i  , j+1, k  )
    i5 = cart2corner(t, i  , j  , k+1)
    i6 = cart2corner(t, i+1, j  , k+1)
    i7 = cart2corner(t, i+1, j+1, k+1)
    i8 = cart2corner(t, i  , j+1, k+1)
    connect((i1, i2, i3, i4, i5, i6, i7, i8), Hexahedron)
  else
    throw(ErrorException("not implemented"))
  end
end

nelements(t::GridTopology) = prod(t.dims)

function facet(t::GridTopology{D}, ind) where {D}
  if D == 1
    ind
  elseif D == 2
    N = 2prod(t.dims)
    if ind ≤ N
      if isodd(ind)
        i, j = elem2cart(t, (ind + 1) ÷ 2)
        i1 = cart2corner(t, i,   j)
        i2 = cart2corner(t, i+1, j)
      else
        i, j = elem2cart(t, ind ÷ 2)
        i1 = cart2corner(t, i,   j)
        i2 = cart2corner(t, i, j+1)
      end
    else
      if isodd(ind)
        i = t.dims[1] + 1
        j = ((ind - N) + 1) ÷ 2
        i1 = cart2corner(t, i,   j)
        i2 = cart2corner(t, i, j+1)
      else
        i = (ind - N) ÷ 2
        j = t.dims[2] + 1
        i1 = cart2corner(t, i,   j)
        i2 = cart2corner(t, i+1, j)
      end
    end
    connect((i1, i2), Segment)
  elseif D == 3
    throw(ErrorException("not implemented"))
  else
    throw(ErrorException("not implemented"))
  end
end

function nfacets(t::GridTopology{D}) where {D}
  if D == 1
    t.dims[1] + 1
  elseif D == 2
    2prod(t.dims) +
    t.dims[1] +
    t.dims[2]
  elseif D == 3
    3prod(t.dims) +
    prod(t.dims[[1,2]]) +
    prod(t.dims[[1,3]]) +
    prod(t.dims[[2,3]])
  else
    throw(ErrorException("not implemented"))
  end
end
