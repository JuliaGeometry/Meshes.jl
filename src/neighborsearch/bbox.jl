# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BVHNode{B}

A node in a bounding volume hierarchy (BVH), storing a bounding box and either
child node indices or a range of positions in the permutation vector.

The fields are defined as follows:

1. `box` stores the bounding box of the node.
2. `a` stores the index of the left child for internal nodes, or the symmetric
   of the starting position in the permutation vector for leaf nodes.
3. `b` stores the index of the right child for internal nodes, or the ending
   position in the permutation vector for leaf nodes.

A negative `a` identifies a leaf node, while a positive `a` identifies an
internal node.
"""
struct BVHNode{B}
  box::B
  a::Int
  b::Int
end

"""
	BBoxSearch(domain, leafsize=8)

A method for searching elements whose bounding boxes intersect the bounding
box of a query using a bounding volume hierarchy (BVH).
"""
struct BBoxSearch{D,B} <: NeighborSearchMethod
  domain::D
  boxes::Vector{B}
  perm::Vector{Int}
  nodes::Vector{BVHNode{B}}
  leafsize::Int
end

function BBoxSearch(domain::D, leafsize=8) where {D<:Domain}
  # validate the leaf size
  leafsize > 0 || throw(ArgumentError("leaf size must be positive"))

  # get the number of elements in the domain and compute their bounding boxes
  n = nelements(domain)

  # compute the bounding boxes of the elements, their corresponding centers, 
  boxes = boundingbox.(domain)
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
  sizehint!(nodes, 2 * cld(n, leafsize))

  # recursively build the BVH and store the nodes in the `nodes` vector
  _buildbvh!(nodes, boxes, centers, perm, 1, n, leafsize)

  # return the constructed BVH
  BBoxSearch{D,B}(domain, boxes, perm, nodes, leafsize)
end

BBoxSearch(geoms, leafsize=8) = BBoxSearch(GeometrySet(geoms), leafsize)

function _buildbvh!(nodes, boxes, centers, perm, first, last, leafsize)
  # reserve the node position so that the root remains node 1.
  nodeind = length(nodes) + 1
  push!(nodes, BVHNode(boxes[perm[first]], 0, 0))

  # arrived at a leaf node, store the range of elements and return; the negative first index indicates a leaf
  if last - first + 1 ≤ leafsize
    nodebox = _bboxes(boxes[perm[i]] for i in first:last)
    nodes[nodeind] = BVHNode(nodebox, -first, last)
    return nodeind
  end

  # get axis to be splitted
  axis = _splitaxis(centers, perm, first, last)

  # sort the elements in the range [first, last] by the center of their bounding boxes along the chosen axis
  range = view(perm, first:last)
  sort!(range; by=i -> centers[i][axis])
  localmiddle = cld(length(range), 2)
  middle = first + localmiddle - 1

  # recursively build and store in `nodes` the left and right children
  left = _buildbvh!(nodes, boxes, centers, perm, first, middle, leafsize)
  right = _buildbvh!(nodes, boxes, centers, perm, middle + 1, last, leafsize)

  # compute the bounding box of the current node by merging the bounding boxes of its left and right children
  nodebox = _bboxes((nodes[left].box, nodes[right].box))

  # store the current node with the bounding box and the indices of the left and right children
  nodes[nodeind] = BVHNode(nodebox, left, right)

  # return the current node index to the caller so that it can be stored in its parent node
  nodeind
end

function _splitaxis(centers, perm, first, last)
  # compute the axis with the widest spread of the centers of the bounding boxes in the range [first, last]
  lo = centers[perm[first]]
  hi = lo

  @inbounds for k in (first + 1):last
    center = centers[perm[k]]

    lo = min.(lo, center)
    hi = max.(hi, center)
  end

  argmax(hi .- lo)
end

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

Apply `f` to the index of each candidate selected by `method` for `geom`.

This function is allocation-free apart from its internal traversal stack and is intended as the primitive 
interface for algorithms that process candidates immediately, avoiding the need to materialize an intermediate
vector of indices.

See also: [`search`](@ref), [`search!`](@ref).
"""
function _foreachcandidate(f, geom, method::BBoxSearch)
  # compute the bounding box of the query geometry `geom` and initialize a stack with the root node index
  geombox = boundingbox(geom)
  stack = [1]

  # traverse the BVH using a stack-based approach to find all nodes whose bounding boxes intersect with the query bounding box
  while !isempty(stack)
    # pop the last node index from the stack and retrieve the corresponding node from the BVH
    nodeid = pop!(stack)
    node = method.nodes[nodeid]

    # check if the bounding box of the current node intersects with the geom bounding box; if not, skip to the next iteration
    intersects(node.box, geombox) || continue

    # if the current node is a leaf node, check each element in the range [first, last] to see if its bounding box intersects
    # with the geom bounding box; if so, add its index to the output vector
    if _isleaf(node)
      @inbounds for k in abs(node.a):node.b
        ind = method.perm[k]
        intersects(method.boxes[ind], geombox) || continue
        f(ind)
      end
    else
      push!(stack, node.a)
      push!(stack, node.b)
    end
  end
end

_isleaf(node::BVHNode) = node.a < 0
