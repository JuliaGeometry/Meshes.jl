# Implementation of Bentley-Ottmann algorith
# https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm

"""
    BentleyOttmann(segments)

Compute pairwise intersections between n `segments`
in O(n⋅log(n)) time using Bentley-Ottmann sweep line
algorithm.

Outputs a Dictionary of {Point, Vector{Tuple{Point, Point}}}
where the key is each intersection point and the values are all
pairs of segments that intersect at that point.
"""
function BentleyOttmann(segments)
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
  lookup = Dict{S,Int}()
  for (i, s) in enumerate(segs)
    a, b = extrema(s)
    BinaryTrees.insert!(𝒬, a)
    BinaryTrees.insert!(𝒬, b)
    haskey(ℒ, a) ? push!(ℒ[a], (a, b)) : (ℒ[a] = [(a, b)])
    haskey(𝒰, b) ? push!(𝒰[b], (a, b)) : (𝒰[b] = [(a, b)])
    haskey(ℒ, b) || (ℒ[b] = S[])
    haskey(𝒰, a) || (𝒰[a] = S[])
    haskey(𝒞, a) || (𝒞[a] = S[])
    haskey(𝒞, b) || (𝒞[b] = S[])
    lookup[(a, b)] = i
  end

  # sweep line
  I = Dict{P,Vector{Tuple{IntersectionType,Vector{Int}}}}()
  while !isnothing(BinaryTrees.root(𝒬))
    p = BinaryTrees.key(BinaryTrees.minnode(𝒬))
    BinaryTrees.delete!(𝒬, p)
    Meshes.handle!(I, lookup, p, S, 𝒬, 𝒯, ℒ, 𝒰, 𝒞)
  end
  I
end

function handle!(I, lookup, p, S, 𝒬, 𝒯, ℒ, 𝒰, 𝒞)
  𝒮ₛ = get(ℒ, p, S[])
  𝒮ₑ = get(𝒰, p, S[])
  𝒮ᵢ = get(𝒞, p, S[])
  _processends!(𝒮ₑ, 𝒬, 𝒯, 𝒞)
  _processstarts!(𝒮ₛ, 𝒬, 𝒯, 𝒞)
  __processintersects!(𝒮ᵢ, 𝒬, 𝒯, 𝒞)
  println("p: ", p)
  println("𝒮ₑ:", 𝒮ₑ)
  println("𝒮ᵢ: ", 𝒮ᵢ)
  println("-----------------")
  if !isempty(𝒮ₛ ∪ 𝒮ₑ ∪ 𝒮ᵢ)
    corners = 𝒮ₛ ∪ 𝒮ₑ
    crossings = 𝒮ᵢ
    I[p] = _intersection(lookup, corners, crossings)
  end
end

function _processstarts!(𝒮ₛ, 𝒬, 𝒯, 𝒞)
  [BinaryTrees.insert!(𝒯, s) for s in 𝒮ₛ]
  for s in 𝒮ₛ
    prev, next = BinaryTrees.prevnext(𝒯, s)
    s = Segment(s)
    if !isnothing(prev) && !isnothing(next)
      new_geom, new_type = _newevent(Segment(BinaryTrees.key(next)), Segment(BinaryTrees.key(prev)))
      if _checkintersection(new_type)
        BinaryTrees.insert!(𝒬, new_geom)
        haskey(𝒞, new_geom) ? push!(𝒞[new_geom], BinaryTrees.key(next), BinaryTrees.key(prev)) :
        (𝒞[new_geom] = [BinaryTrees.key(next), BinaryTrees.key(prev)])
      end
    end
    if !isnothing(prev)
      new_geom, new_type = _newevent(Segment(BinaryTrees.key(prev)), s)
      if new_type == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
        haskey(𝒞, new_geom) ? push!(𝒞[new_geom], BinaryTrees.key(prev), vertices(s)) :
        (𝒞[new_geom] = [BinaryTrees.key(prev), vertices(s)])
      end
    end
    if !isnothing(next)
      new_geom, new_type = _newevent(s, Segment(BinaryTrees.key(next)))
      if new_type == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
        haskey(𝒞, new_geom) ? push!(𝒞[new_geom], vertices(s), BinaryTrees.key(next)) :
        (𝒞[new_geom] = [vertices(s), BinaryTrees.key(next)])
      end
    end
  end
end

function _processends!(𝒮ₑ, 𝒬, 𝒯, 𝒞)
  for s in 𝒮ₑ
    prev, next = BinaryTrees.prevnext(𝒯, s)
    BinaryTrees.delete!(𝒯, s)
    s = Segment(s)
    if !isnothing(prev) && !isnothing(next)
      new_geom, new_type = _newevent(Segment(BinaryTrees.key(next)), Segment(BinaryTrees.key(prev)))
      if new_type == IntersectionType(0)
        BinaryTrees.insert!(𝒬, new_geom)
        haskey(𝒞, new_geom) ? push!(𝒞[new_geom], BinaryTrees.key(next), BinaryTrees.key(prev)) :
        (𝒞[new_geom] = [BinaryTrees.key(next), BinaryTrees.key(prev)])
      end
    end
  end
end

function __processintersects!(𝒮ᵢ, 𝒬, 𝒯, 𝒞)
  for s in 𝒮ᵢ
    prev, _ = BinaryTrees.prevnext(𝒯, s)
    if !isnothing(prev)

      # Find segments r and u
      r, _ = BinaryTrees.prevnext(𝒯, s)
      _, u = BinaryTrees.prevnext(𝒯, BinaryTrees.key(prev))

      # Remove crossing points rs and tu from event queue
      if !isnothing(r)
        new_geom, new_type = _newevent(Segment(BinaryTrees.key(r)), Segment(s))
        if new_type == IntersectionType(0)
          BinaryTrees.delete!(𝒬, new_geom)
        end
      end
      if !isnothing(u)
        new_geom, new_type = _newevent(Segment(BinaryTrees.key(u)), Segment(BinaryTrees.key(prev)))
        if new_type == IntersectionType(0)
          BinaryTrees.delete!(𝒬, new_geom)
        end
      end

      # Add crossing points rt and su to event queue
      if !isnothing(r)
        new_geom, new_type = _newevent(Segment(BinaryTrees.key(r)), Segment(BinaryTrees.key(prev)))
        if new_type == IntersectionType(0)
          BinaryTrees.insert!(𝒬, new_geom)
          haskey(𝒞, new_geom) ? push!(𝒞[new_geom], BinaryTrees.key(r), BinaryTrees.key(prev)) :
          (𝒞[new_geom] = [BinaryTrees.key(r), BinaryTrees.key(prev)])
        end
      end
      if !isnothing(u)
        new_geom, new_type = _newevent(Segment(BinaryTrees.key(u)), Segment(s))
        if new_type == IntersectionType(0)
          BinaryTrees.insert!(𝒬, new_geom)
          haskey(𝒞, new_geom) ? push!(𝒞[new_geom], BinaryTrees.key(u), s) : (𝒞[new_geom] = [BinaryTrees.key(u), s])
        end
      end
    end
  end
end

function _pushintersection(lookup, corners, crossings)
  return [
    (IntersectionType(4), [lookup[segment] for segment in corners]),
    (IntersectionType(0), [lookup[segment] for segment in crossings])
  ]
end
_key(node::Nothing) = nothing
function _newevent(s₁::Segment, s₂::Segment)
  new_event = intersection(s₁, s₂)
  if !isnothing(new_event)
    get(new_event), type(new_event)
  else
    nothing, nothing
  end
end

function _checkintersection(type)
  type == IntersectionType(0) || type == IntersectionType(1) || type == IntersectionType(2)
end
