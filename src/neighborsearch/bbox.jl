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

BBoxSearch(boxes::Vector{<:Box}, leafsize=8) = BBoxSearch(BVHTree(copy(boxes), leafsize))

function BBoxSearch(domain::Domain, leafsize=8)
  boxes = boundingbox.(domain)
  BBoxSearch(BVHTree(boxes, leafsize))
end

BBoxSearch(geoms, leafsize=8) = BBoxSearch(GeometrySet(geoms), leafsize)

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
