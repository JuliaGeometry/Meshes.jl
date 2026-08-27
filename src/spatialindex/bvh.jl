# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BVH(domain; leafsize=8)

Construct a static bounding volume hierarchy over the elements of
`domain`.

The hierarchy stores axis-aligned bounding boxes and supports broad-phase
queries with [`candidates`](@ref) and [`candidates!`](@ref).

The `leafsize` parameter specifies the maximum number of elements stored
in each leaf node. The root node is always the first node in the `nodes` 
vector.
"""
struct BVHNode{B}
  box::B
  left::Int
  right::Int
  first::Int
  last::Int
end

struct BVH{D,B} <: SpatialIndex
  domain::D
  boxes::Vector{B}
  perm::Vector{Int}
  nodes::Vector{BVHNode{B}}
  leafsize::Int
end

function BVH(domain::Domain; leafsize::Int=8)
  # validate the leaf size
  leafsize > 0 || throw(ArgumentError("leaf size must be positive"))

  # get the number of elements in the domain and compute their bounding boxes
  n = nelements(domain)

  # compute the bounding boxes of the elements, their corresponding centers, 
  boxes = [boundingbox(element(domain, i)) for i in 1:n]
  centers = map(boxes) do box
    c = coords(center(box))
    names = propertynames(c)
    ntuple(i -> getproperty(c, names[i]), length(names))
  end

  # initialize a permutation vector
  perm = collect(1:n)

  # determine the element type of the bounding boxes and initialize an empty vector of BVH nodes
  B = eltype(boxes)
  nodes = BVHNode{B}[]

  # recursively build the BVH and store the nodes in the `nodes` vector
  _buildbvh!(nodes, boxes, centers, perm, 1, n, leafsize)

  # return the constructed BVH
  BVH{typeof(domain),B}(domain, boxes, perm, nodes, leafsize)
end

function _buildbvh!(nodes, boxes, centers, perm, first, last, leafsize)
  # reserve the node position so that the root remains node 1.
  nodeind = length(nodes) + 1
  push!(nodes, BVHNode(boxes[perm[first]], 0, 0, 0, 0))

  # arrived at a leaf node, store the range of elements and return
  if last - first + 1 ≤ leafsize
    nodebox = _bboxes(boxes[perm[i]] for i in first:last)
    nodes[nodeind] = BVHNode(nodebox, 0, 0, first, last)
    return nodeind
  end

  # get axis to be splitted
  axis = _splitaxis(centers, perm, first, last)

  # partial sort the elements in the range [first, last] by the center of their bounding boxes along the chosen axis
  range = view(perm, first:last)
  localmiddle = cld(length(range), 2)
  partialsort!(range, localmiddle; by=i -> centers[i][axis]) # range is a view of perm, so perm is modified in place!
  middle = first + localmiddle - 1

  # recursively build and store in `nodes` the left and right children
  left = _buildbvh!(nodes, boxes, centers, perm, first, middle, leafsize)
  right = _buildbvh!(nodes, boxes, centers, perm, middle + 1, last, leafsize)

  # compute the bounding box of the current node by merging the bounding boxes of its left and right children
  nodebox = _bboxes((nodes[left].box, nodes[right].box))

  # store the current node with the bounding box and the indices of the left and right children
  # the first and last indices are not used for non-leaf nodes, so they are set to 0
  nodes[nodeind] = BVHNode(nodebox, left, right, 0, 0)

  # return the current node index to the caller so that it can be stored in its parent node
  nodeind
end

function _splitaxis(centers, perm, first, last)
  # compute the axis with the widest spread of the centers of the bounding boxes in the range [first, last]
  dim = length(centers[perm[first]])

  lo = collect(centers[perm[first]])
  hi = copy(lo)

  for k in (first + 1):last
    center = centers[perm[k]]
    for axis in 1:dim
      value = center[axis]
      lo[axis] = min(lo[axis], value)
      hi[axis] = max(hi[axis], value)
    end
  end

  return argmax(hi .- lo)
end

function candidates(query, bvh::BVH)
  inds = Int[]
  candidates!(inds, query, bvh)
end

"""
    candidates!(inds, query, bvh)
Find the indices of the elements in `bvh` whose bounding boxes intersect with the bounding box of `query`. 
The indices are stored in the preallocated vector `inds`, which is cleared before storing the results. 
The function returns the vector `inds` for convenience.
"""
function candidates!(inds::Vector{Int}, query, bvh::BVH)
  # clear the output vector to ensure it only contains the results of the current query
  empty!(inds)

  # compute the bounding box of the query and initialize a stack with the root node index
  querybox = boundingbox(query)
  stack = Int[1]

  # traverse the BVH using a stack-based approach to find all nodes whose bounding boxes intersect with the query bounding box
  while !isempty(stack)
    # pop the last node index from the stack and retrieve the corresponding node from the BVH
    nodeid = pop!(stack)
    node = bvh.nodes[nodeid]

    # check if the bounding box of the current node intersects with the query bounding box; if not, skip to the next iteration
    intersects(node.box, querybox) || continue

    # if the current node is a leaf node, check each element in the range [first, last] to see if its bounding box intersects
    # with the query bounding box; if so, add its index to the output vector
    if _isleaf(node)
      for k in node.first:node.last
        ind = bvh.perm[k]
        if intersects(bvh.boxes[ind], querybox)
          push!(inds, ind)
        end
      end
    else
      push!(stack, node.left)
      push!(stack, node.right)
    end
  end

  return inds
end

_isleaf(node::BVHNode) = iszero(node.left)
