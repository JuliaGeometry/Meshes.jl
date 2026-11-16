# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BVHNode{B}

A node in a bounding volume hierarchy (BVH), storing a bounding box and either
child node indices or a range of positions in the permutation vector.

The fields are defined as follows:

1. `box` stores the bounding box of the node.
2. `left` stores the index of the left child for internal nodes, or the symmetric
   of the starting position in the permutation vector for leaf nodes.
3. `right` stores the index of the right child for internal nodes, or the ending
   position in the permutation vector for leaf nodes.

A negative `left` identifies a leaf node, while a positive `left` identifies an
internal node.
"""
struct BVHNode{B<:Box}
  box::B
  left::Int
  right::Int
end

"""
    isleaf(node)

Return `true` if `node` is a leaf node, `false` otherwise.
"""
isleaf(node::BVHNode) = node.left < 0

"""
    BVHTree(boxes; leafsize=8)

Construct a static bounding volume hierarchy over the bounding boxes in `boxes`.
The positive integer `leafsize` specifies the maximum number of boxes per leaf.

The input vector is stored without copying and must not be modified while
the tree is in use.
"""
struct BVHTree{B<:Box}
  boxes::Vector{B}
  perm::Vector{Int}
  nodes::Vector{BVHNode{B}}
  leafsize::Int
end

function BVHTree(boxes::AbstractVector{<:Box}; leafsize=8)
  # sanity check
  leafsize > 0 || throw(ArgumentError("leaf size must be positive"))

  # number of boxes
  nboxes = length(boxes)

  # initialize an empty vector of BVH nodes
  nodes = BVHNode{eltype(boxes)}[]
  sizehint!(nodes, _bvhnnodes(nboxes, leafsize))

  # compute the centers of the bounding boxes
  centers = map(boxes) do box
    CoordRefSystems.values(coords(center(box)))
  end

  # map nodes in the tree to original box indices
  perm = collect(1:nboxes)

  # recursively build the BVH and store the nodes in the nodes vector
  _bvhbuild!(nodes, boxes, centers, perm, 1, nboxes, leafsize)

  # return the constructed BVH
  BVHTree(boxes, perm, nodes, leafsize)
end

"""
    foreachcandidate(func, box, tree, nodeind=1)

Apply `func` to each original element index whose cached box in `tree` intersects `box`.

This function is allocation-free and is intended as the primitive interface for algorithms that 
process candidates immediately, avoiding the need to materialize an intermediate vector of indices.
"""
function foreachcandidate(func, box::Box, tree::BVHTree, nodeind=1)
  node = tree.nodes[nodeind]

  # early exit if the current node's box does not intersect with the query box
  intersects(node.box, box) || return nothing

  # if the current node is a leaf node, check each element in the range
  # [firstind, lastind] to see if its box intersects with the query box;
  # if so, apply the function
  if isleaf(node)
    @inbounds for k in abs(node.left):node.right
      ind = tree.perm[k]
      intersects(tree.boxes[ind], box) || continue
      func(ind)
    end
  else
    foreachcandidate(func, box, tree, node.right)
    foreachcandidate(func, box, tree, node.left)
  end

  nothing
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _bvhnnodes(nboxes, leafsize)
  width = prevpow(2, cld(nboxes, leafsize))
  n, p = divrem(nboxes, width)
  splits = if n > leafsize
    width
  elseif n == leafsize
    p
  else
    0
  end
  2 * (width + splits) - 1
end

function _bvhbuild!(nodes, boxes, centers, perm, firstind, lastind, leafsize)
  # reserve the current node before constructing its descendants
  ind = length(nodes) + 1
  push!(nodes, BVHNode(boxes[perm[firstind]], 0, 0))
  leafnode = lastind - firstind + 1 ≤ leafsize

  if leafnode
    # this is a leaf node, compute the bounding box for the range and store the leaf range
    box = boundingbox(boxes[perm[i]] for i in firstind:lastind)
    left = -firstind
    right = lastind
  else
    # this is an internal node, split and recursively build the children
    dim = _bvhdim(centers, perm, firstind, lastind)

    # sort the elements in the range [firstind, lastind] by the center of
    # their bounding boxes along the split dimension
    pview = view(perm, firstind:lastind)
    sort!(pview; by=i -> centers[i][dim])
    localind = cld(length(pview), 2)
    middleind = firstind + localind - 1

    # recursively build and store in `nodes` the left and right children
    left = _bvhbuild!(nodes, boxes, centers, perm, firstind, middleind, leafsize)
    right = _bvhbuild!(nodes, boxes, centers, perm, middleind + 1, lastind, leafsize)

    # current nodes' bounding box is the bounding box of its left and right children
    box = boundingbox((nodes[left].box, nodes[right].box))
  end

  # store the current node with its bounding box and either the left and right children indices or the leaf range
  nodes[ind] = BVHNode(box, left, right)

  # return the current node index to the caller so that it can be stored in its parent node
  ind
end

function _bvhdim(centers, perm, firstind, lastind)
  # dimension with the widest spread of the centers of
  # the bounding boxes in the range [firstind, lastind]
  lo = centers[perm[firstind]]
  hi = lo
  @inbounds for k in (firstind + 1):lastind
    cₖ = centers[perm[k]]
    lo = min.(lo, cₖ)
    hi = max.(hi, cₖ)
  end
  argmax(hi .- lo)
end
