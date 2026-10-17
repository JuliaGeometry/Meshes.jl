# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BBoxSearch(domain, leafsize=8)

A method for searching elements whose bounding boxes intersect the bounding
box of a query using a bounding volume hierarchy (BVH).
"""
struct BBoxSearch{T<:BVHTree} <: NeighborSearchMethod
  tree::T
end

function BBoxSearch(domain::Domain, leafsize=8)
  BBoxSearch(BVHTree(domain, leafsize))
end

BBoxSearch(geoms, leafsize=8) = BBoxSearch(GeometrySet(geoms), leafsize)

function search(geom, method::BBoxSearch; mask=nothing)
  inds = Int[]
  search!(inds, geom, method; mask=mask)
end

function search!(neighbors, geom, method::BBoxSearch; mask=nothing)
  # clear the output vector to ensure it only contains the results of the current query geometry
  empty!(neighbors)

  _foreachcandidate(geom, method) do ind
    if isnothing(mask) || mask[ind]
      push!(neighbors, ind)
    end
  end

  neighbors
end

"""
    _foreachcandidate(f, geom, method)

Apply `f` to each original candidate index for `geom` using the tree of `method`.
"""
_foreachcandidate(f, geom, method::BBoxSearch) = _foreachcandidate(f, boundingbox(geom), method.tree)
