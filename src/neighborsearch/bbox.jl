# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BBoxSearch(domain; leafsize=8)

A method for searching elements whose bounding boxes intersect the bounding
box of a query using a bounding volume hierarchy (BVH).
"""
struct BBoxSearch{T<:BVHTree} <: NeighborSearchMethod
  tree::T
end

function BBoxSearch(domain::Domain; leafsize=8)
  boxes = boundingbox.(domain)
  tree = BVHTree(boxes; leafsize)
  BBoxSearch(tree)
end

BBoxSearch(geoms; leafsize=8) = BBoxSearch(GeometrySet(geoms); leafsize)

function search(geom, method::BBoxSearch; mask=nothing)
  inds = Int[]
  stack = Int[]

  foreachcandidate(boundingbox(geom), method.tree, stack) do ind
    if isnothing(mask) || mask[ind]
      push!(inds, ind)
    end
  end

  inds
end
