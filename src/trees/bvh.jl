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
  sizehint!(nodes, 2 * cld(nboxes, leafsize))

  # compute the centers of the bounding boxes
  centers = map(boxes) do box
    CoordRefSystems.values(coords(center(box)))
  end

  # map nodes in the tree to original box indices
  perm = collect(1:nboxes)

  # recursively build the BVH and store the nodes in the nodes vector
  _buildbvh!(nodes, boxes, centers, perm, 1, nboxes, leafsize)

  # return the constructed BVH
  BVHTree(boxes, perm, nodes, leafsize)
end

function _buildbvh!(nodes, boxes, centers, perm, firstind, lastind, leafsize)
  # reserve the node position so that the root remains node 1.
  ind = length(nodes) + 1
  box = boxes[perm[firstind]]
  push!(nodes, BVHNode(box, 0, 0))

  # arrived at a leaf node, store the range of elements and return; the negative first index indicates a leaf
  if lastind - firstind + 1 ≤ leafsize
    box = boundingbox(boxes[perm[i]] for i in firstind:lastind)
    nodes[ind] = BVHNode(box, -firstind, lastind)
    return ind
  end

  # find dimension to split
  dim = _splitdim(centers, perm, firstind, lastind)

  # sort the elements in the range [firstind, lastind] by the center of their bounding boxes along the chosen dim
  pview = view(perm, firstind:lastind)
  sort!(pview; by=i -> centers[i][dim])
  localind = cld(length(pview), 2)
  middleind = firstind + localind - 1

  # recursively build and store in `nodes` the left and right children
  left = _buildbvh!(nodes, boxes, centers, perm, firstind, middleind, leafsize)
  right = _buildbvh!(nodes, boxes, centers, perm, middleind + 1, lastind, leafsize)

  # store the current node with the bounding box and the indices of the left and right children
  box = boundingbox((nodes[left].box, nodes[right].box))
  nodes[ind] = BVHNode(box, left, right)

  # return the current node index to the caller so that it can be stored in its parent node
  ind
end

function _splitdim(centers, perm, firstind, lastind)
  # compute the axis with the widest spread of the centers of the bounding boxes in the range [firstind, lastind]
  lo = centers[perm[firstind]]
  hi = lo

  @inbounds for k in (firstind + 1):lastind
    center = centers[perm[k]]

    lo = min.(lo, center)
    hi = max.(hi, center)
  end

  argmax(hi .- lo)
end

"""
    foreachcandidate(f, geombox, tree, stack)

Apply `f` to each original element index whose cached box in `tree` intersects `geombox`.

This function is allocation-free apart from its internal traversal stack and is intended as the primitive 
interface for algorithms that process candidates immediately, avoiding the need to materialize an intermediate
vector of indices.
"""
function foreachcandidate(f::F, geombox::Box, tree::BVHTree, stack) where {F}
  # initialize a stack with the root node index
  empty!(stack)
  push!(stack, 1)

  # traverse the BVH using a stack-based approach to find all nodes whose bounding boxes intersect with the query bounding box
  while !isempty(stack)
    # pop the last node index from the stack and retrieve the corresponding node from the BVH
    nodeid = pop!(stack)
    node = tree.nodes[nodeid]

    # check if the bounding box of the current node intersects with the geom bounding box; if not, skip to the next iteration
    intersects(node.box, geombox) || continue

    # if the current node is a leaf node, check each element in the range [firstind, lastind] to see if its bounding box intersects
    # with the geom bounding box; if so, add its index to the output vector
    if isleaf(node)
      @inbounds for k in abs(node.left):node.right
        ind = tree.perm[k]
        intersects(tree.boxes[ind], geombox) || continue
        f(ind)
      end
    else
      push!(stack, node.left)
      push!(stack, node.right)
    end
  end
end

isleaf(node::BVHNode) = node.left < 0
