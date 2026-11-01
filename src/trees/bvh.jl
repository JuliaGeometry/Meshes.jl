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
struct BVHNode{B}
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
struct BVHTree{B}
  boxes::Vector{B}
  perm::Vector{Int}
  nodes::Vector{BVHNode{B}}
  leafsize::Int
end

function BVHTree(boxes::AbstractVector{<:Box}; leafsize=8)
  # validate the leaf size
  leafsize > 0 || throw(ArgumentError("leaf size must be positive"))

  n = length(boxes)

  # compute the centers of the bounding boxes
  centers = map(boxes) do box
    CoordRefSystems.values(coords(center(box)))
  end

  # initialize a permutation vector to map positions in the tree partition to original box indices
  perm = collect(1:n)

  # determine the element type of the bounding boxes and initialize an empty vector of BVH nodes
  B = eltype(boxes)
  nodes = BVHNode{B}[]
  sizehint!(nodes, 2 * cld(n, leafsize))

  # recursively build the BVH and store the nodes in the `nodes` vector
  _buildbvh!(nodes, boxes, centers, perm, 1, n, leafsize)

  # return the constructed BVH
  BVHTree(boxes, perm, nodes, leafsize)
end

function _buildbvh!(nodes, boxes, centers, perm, firstind, lastind, leafsize)
  # reserve the node position so that the root remains node 1.
  nodeind = length(nodes) + 1
  push!(nodes, BVHNode(boxes[perm[firstind]], 0, 0))

  # arrived at a leaf node, store the range of elements and return; the negative first index indicates a leaf
  if lastind - firstind + 1 ≤ leafsize
    nodebox = boundingbox(boxes[perm[i]] for i in firstind:lastind)
    nodes[nodeind] = BVHNode(nodebox, -firstind, lastind)
    return nodeind
  end

  # get axis to be splitted
  axis = _splitaxis(centers, perm, firstind, lastind)

  # sort the elements in the range [firstind, lastind] by the center of their bounding boxes along the chosen axis
  range = view(perm, firstind:lastind)
  sort!(range; by=i -> centers[i][axis])
  localmiddle = cld(length(range), 2)
  middle = firstind + localmiddle - 1

  # recursively build and store in `nodes` the left and right children
  left = _buildbvh!(nodes, boxes, centers, perm, firstind, middle, leafsize)
  right = _buildbvh!(nodes, boxes, centers, perm, middle + 1, lastind, leafsize)

  # compute the bounding box of the current node by merging the bounding boxes of its left and right children
  nodebox = boundingbox((nodes[left].box, nodes[right].box))

  # store the current node with the bounding box and the indices of the left and right children
  nodes[nodeind] = BVHNode(nodebox, left, right)

  # return the current node index to the caller so that it can be stored in its parent node
  nodeind
end

function _splitaxis(centers, perm, firstind, lastind)
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
