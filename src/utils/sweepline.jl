# Implementation of Bentley-Ottmann algorith
# https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm

using BinaryTrees


"""
    bentleyottmann(segments)

Compute pairwise intersections between n `segments`
in O(n⋅log(n)) time using Bentley-Ottmann sweep line
algorithm.
"""
function bentleyottmann(segments)
  # adjust vertices of segments
  segs = map(segments) do s
    a, b = extrema(s)
    a > b ? reverse(s) : s
  end

  # retrieve relevant info
  s = first(segs)
  p = minimum(s)
  P = typeof(p)
  S = Tuple{P,P}

  # initialization
  𝒬 = BinaryTrees.AVLTree{P}()
  𝒯 = BinaryTrees.AVLTree{S}()
  ℒ = Dict{P,Vector{S}}()
  𝒰 = Dict{P,Vector{S}}()
  𝒞 = Dict{P,Vector{S}}()
  for s in segs
    a, b = extrema(s)
    BinaryTrees.insert!(𝒬, a)
    BinaryTrees.insert!(𝒬, b)
    haskey(ℒ, a) ? push!(ℒ[a], (a, b)) : (ℒ[a] = [(a, b)])
    haskey(𝒰, b) ? push!(𝒰[b], (a, b)) : (𝒰[b] = [(a, b)])
    haskey(ℒ, b) || (ℒ[b] = S[])
    haskey(𝒰, a) || (𝒰[a] = S[])
    haskey(𝒞, a) || (𝒞[a] = S[])
    haskey(𝒞, b) || (𝒞[b] = S[])
  end
  m = Point(-Inf, -Inf)
  M = Point(Inf, Inf)
  BinaryTrees.insert!(𝒯, (m, m))
  BinaryTrees.insert!(𝒯, (M, M))

  # sweep line
  I = Dict{P,Vector{S}}()
  while !isnothing(BinaryTrees.root(𝒬))
    p = _key(BinaryTrees.root(𝒬))
    BinaryTrees.delete!(𝒬, p)
    handle!(I, p, 𝒬, 𝒯, ℒ, 𝒰, 𝒞)
  end
  I
end

function handle!(I, p, 𝒬, 𝒯, ℒ, 𝒰, 𝒞)
  # Segments that start, end, or intersect at p
  start_segments = ℒ[p]
  end_segments = 𝒰[p]
  intersection_segments = 𝒞[p]

  # If there are multiple segments intersecting at p, record the intersection
  if length(start_segments ∪ end_segments ∪ intersection_segments) > 1
    I[p] = start_segments ∪ end_segments ∪ intersection_segments
  end

  # Remove segments that end at p from the status structure
  for s in end_segments ∪ intersection_segments
    BinaryTrees.delete!(𝒯, s)
  end

  # Insert segments that start at p into the status structure
  for s in start_segments ∪ intersection_segments
    BinaryTrees.insert!(𝒯, s)
  end
  node = BinaryTrees.root(𝒬)

  # Find new event points caused by the insertion or deletion of segments
  for s in start_segments
    s = Segment(s)
    pred = BinaryTrees.left(node)
    succ = BinaryTrees.right(node)
    ns = Segment(pred, succ)
    if pred !== nothing
      new_geom, new_type = _newevent(s, ns)
      if new_geom == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
      end
    end
    if succ !== nothing
      new_geom, new_type = _newevent(s, ns)
      if new_type == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
      end
    end
  end

  for s in end_segments
    s = Segment(s)
    pred = BinaryTrees.left(node)
    succ = BinaryTrees.right(node)
    ns = Segment(pred, succ)
    if pred !== nothing && succ !== nothing
      new_geom, new_type = _newevent(s, ns)
      if new_type == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
      end
    end
  end
end

_key(node::BinaryTrees.AVLNode) = node.key
_geom(intersect::Intersection) = intersect.geom
_type(intersect::Intersection) = intersect.type
function _newevent(s₁::Segment, s₂::Segment)
  new_event = intersection(s₁, s₂)
  _geom(new_event), _type(new_event)
end


function Segment(node₁::BinaryTrees.BinaryNode, node₂::BinaryTrees.BinaryNode)
  node₁ = _key(node₁)
  node₂ = _key(node₂)
  Segment((node₁, node₂))
end
