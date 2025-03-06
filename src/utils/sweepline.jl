# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
    bentleyottmann(segments)

Compute pairwise intersections between n `segments`
in O(n⋅log(n)) time using Bentley-Ottmann sweep line
algorithm.

Outputs a Dictionary of {Point, Vector{Tuple{Point, Point}}}
where the key is each intersection point and the values are all
pairs of segments that intersect at that point.

## References

* Bentley, J. L., & Ottmann, T. 1979. [Algorithms for reporting and counting geometric intersections]
  (https://www.itseng.org/research/papers/topics/VLSI_Physical_Design_Automation/Physical_Verification/DRC/Geometric_Intersection_Problems/1979-Bentley.pdf)
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
  ℬ = get(ℒ, p, S[])
  ℰ = get(𝒰, p, S[])
  ℐ = get(𝒞, p, S[])
  _processend!(ℰ, 𝒬, 𝒯, 𝒞)
  _processbegin!(ℬ, 𝒬, 𝒯, 𝒞)
  _processintersects!(ℐ, 𝒬, 𝒯, 𝒞)
  if !isempty(ℬ ∪ ℰ ∪ ℐ)
    corners = ℬ ∪ ℰ
    crossings = ℐ
    I[p] = _pushintersection(lookup, corners, crossings)
  end
end

function _processbegin!(ℬ, 𝒬, 𝒯, 𝒞)
  [BinaryTrees.insert!(𝒯, s) for s in ℬ]
  for s in ℬ
    prev, next = BinaryTrees.prevnext(𝒯, s)
    s = Segment(s)
    if !isnothing(prev) && !isnothing(next)
      newgeom, newtype = _newevent(Segment(BinaryTrees.key(next)), Segment(BinaryTrees.key(prev)))
      if _checkintersection(newtype)
        BinaryTrees.insert!(𝒬, newgeom)
        haskey(𝒞, newgeom) ? push!(𝒞[newgeom], BinaryTrees.key(next), BinaryTrees.key(prev)) :
        (𝒞[newgeom] = [BinaryTrees.key(next), BinaryTrees.key(prev)])
      end
    end
    if !isnothing(prev)
      newgeom, newtype = _newevent(Segment(BinaryTrees.key(prev)), s)
      if _checkintersection(newtype)
        BinaryTrees.insert!(𝒬, newgeom)
        haskey(𝒞, newgeom) ? push!(𝒞[newgeom], BinaryTrees.key(prev), vertices(s)) :
        (𝒞[newgeom] = [BinaryTrees.key(prev), vertices(s)])
      end
    end
    if !isnothing(next)
      newgeom, newtype = _newevent(s, Segment(BinaryTrees.key(next)))
      if _checkintersection(newtype)
        BinaryTrees.insert!(𝒬, newgeom)
        haskey(𝒞, newgeom) ? push!(𝒞[newgeom], vertices(s), BinaryTrees.key(next)) :
        (𝒞[newgeom] = [vertices(s), BinaryTrees.key(next)])
      end
    end
  end
end

function _processend!(ℰ, 𝒬, 𝒯, 𝒞)
  for s in ℰ
    prev, next = BinaryTrees.prevnext(𝒯, s)
    BinaryTrees.delete!(𝒯, s)
    s = Segment(s)
    if !isnothing(prev) && !isnothing(next)
      newgeom, newtype = _newevent(Segment(BinaryTrees.key(next)), Segment(BinaryTrees.key(prev)))
      if _checkintersection(newtype)
        BinaryTrees.insert!(𝒬, newgeom)
        haskey(𝒞, newgeom) ? push!(𝒞[newgeom], BinaryTrees.key(next), BinaryTrees.key(prev)) :
        (𝒞[newgeom] = [BinaryTrees.key(next), BinaryTrees.key(prev)])
      end
    end
  end
end

function _processintersects!(ℐ, 𝒬, 𝒯, 𝒞)
  for s in ℐ
    prev, _ = BinaryTrees.prevnext(𝒯, s)
    if !isnothing(prev)

      # Find segments r and u
      r, _ = BinaryTrees.prevnext(𝒯, s)
      _, u = BinaryTrees.prevnext(𝒯, BinaryTrees.key(prev))

      # Remove crossing points rs and tu from event queue
      if !isnothing(r)
        newgeom, newtype = _newevent(Segment(BinaryTrees.key(r)), Segment(s))
        if _checkintersection(newtype)
          BinaryTrees.delete!(𝒬, newgeom)
        end
      end
      if !isnothing(u)
        newgeom, newtype = _newevent(Segment(BinaryTrees.key(u)), Segment(BinaryTrees.key(prev)))
        if _checkintersection(newtype)
          BinaryTrees.delete!(𝒬, newgeom)
        end
      end

      # Add crossing points rt and su to event queue
      if !isnothing(r)
        newgeom, newtype = _newevent(Segment(BinaryTrees.key(r)), Segment(BinaryTrees.key(prev)))
        if _checkintersection(newtype)
          BinaryTrees.insert!(𝒬, newgeom)
          haskey(𝒞, newgeom) ? push!(𝒞[newgeom], BinaryTrees.key(r), BinaryTrees.key(prev)) :
          (𝒞[newgeom] = [BinaryTrees.key(r), BinaryTrees.key(prev)])
        end
      end
      if !isnothing(u)
        newgeom, newtype = _newevent(Segment(BinaryTrees.key(u)), Segment(s))
        if _checkintersection(newtype)
          BinaryTrees.insert!(𝒬, newgeom)
          haskey(𝒞, newgeom) ? push!(𝒞[newgeom], BinaryTrees.key(u), s) : (𝒞[newgeom] = [BinaryTrees.key(u), s])
        end
      end
    end
  end
end

function _pushintersection(lookup, corners, crossings)
  [
    (CornerCrossing, unique([lookup[segment] for segment in corners])),
    (Crossing, unique([lookup[segment] for segment in crossings]))
  ]
end
function _newevent(s₁::Segment, s₂::Segment)
  newevent = intersection(s₁, s₂)
  if !isnothing(newevent)
    get(newevent), type(newevent)
  else
    nothing, nothing
  end
end

function _checkintersection(type)
  type == Crossing || type == EdgeTouching
end
