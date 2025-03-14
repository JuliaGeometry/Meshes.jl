# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    bentleyottmann(segments)

Compute pairwise intersections between n `segments` in
O(n⋅log(n)) time using Bentley-Ottmann sweep line algorithm.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function bentleyottmann(segments)
  # orient segments
  segs = map(segments) do s
    a, b = extrema(s)
    a > b ? (b, a) : (a, b)
  end

  # retrieve types
  P = eltype(first(segs))
  S = Tuple{P,P}

  # initialization
  𝒬 = BinaryTrees.AVLTree{P}()
  ℛ = BinaryTrees.AVLTree{S}()
  ℬ = Dict{P,Vector{S}}()
  ℰ = Dict{P,Vector{S}}()
  ℳ = Dict{P,Vector{S}}()
  lookup = Dict{S,Int}()
  for (i, (a, b)) in enumerate(segs)
    BinaryTrees.insert!(𝒬, a)
    BinaryTrees.insert!(𝒬, b)
    haskey(ℬ, a) ? push!(ℬ[a], (a, b)) : (ℬ[a] = [(a, b)])
    haskey(ℰ, b) ? push!(ℰ[b], (a, b)) : (ℰ[b] = [(a, b)])
    lookup[(a, b)] = i
  end

  # sweep line algorithm
  points = Vector{P}()
  seginds = Vector{Vector{Int}}()
  while !_isempty(𝒬)
    p = BinaryTrees.key(BinaryTrees.minnode(𝒬))
    BinaryTrees.delete!(𝒬, p)
    _handle!(points, seginds, lookup, p, 𝒬, ℛ, ℬ, ℰ, ℳ)
  end
  points, seginds
end

function _handle!(points, seginds, lookup, p, 𝒬, ℛ, ℬ, ℰ, ℳ)
  ℬₚ = _segmentswith(ℬ, p)
  ℰₚ = _segmentswith(ℰ, p)
  ℳₚ = _segmentswith(ℳ, p)
  _processend!(ℰₚ, 𝒬, ℛ, ℳ)
  _processbeg!(ℬₚ, 𝒬, ℛ, ℳ)
  _processmid!(ℳₚ, 𝒬, ℛ, ℳ)
  inds = [lookup[s] for s in ℬₚ ∪ ℰₚ ∪ ℳₚ]
  if !isempty(inds)
    push!(points, p)
    push!(seginds, inds)
  end
end

function _processbeg!(ℬₚ, 𝒬, ℛ, ℳ)
  for s in ℬₚ
    BinaryTrees.insert!(ℛ, s)
  end
  for s in ℬₚ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev)
      _newevent!(𝒬, ℳ, BinaryTrees.key(prev), s)
    end
    if !isnothing(next)
      _newevent!(𝒬, ℳ, s, BinaryTrees.key(next))
    end
  end
end

function _processend!(ℰₚ, 𝒬, ℛ, ℳ)
  for s in ℰₚ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev) && !isnothing(next)
      _newevent!(𝒬, ℳ, BinaryTrees.key(next), BinaryTrees.key(prev))
    end
    BinaryTrees.delete!(ℛ, s)
  end
end

function _processmid!(ℳₚ, 𝒬, ℛ, ℳ)
  for s in ℳₚ
    prev, _ = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev)
      # find segments r and u
      r = prev
      _, u = BinaryTrees.prevnext(ℛ, BinaryTrees.key(prev))

      # remove crossing points rs and tu from event queue
      if !isnothing(r)
        _rmevent!(𝒬, BinaryTrees.key(r), s)
      end
      if !isnothing(u)
        _rmevent!(𝒬, BinaryTrees.key(u), BinaryTrees.key(prev))
      end

      # add crossing points rt and su to event queue
      if !isnothing(r)
        _newevent!(𝒬, ℳ, BinaryTrees.key(r), BinaryTrees.key(prev))
      end
      if !isnothing(u)
        _newevent!(𝒬, ℳ, BinaryTrees.key(u), s)
      end
    end
  end
end

function _newevent!(𝒬, ℳ, s₁, s₂)
  intersection(Segment(s₁), Segment(s₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      p = get(I)
      BinaryTrees.insert!(𝒬, p)
      if haskey(ℳ, p)
        push!(ℳ[p], s₁, s₂)
      else
        ℳ[p] = [s₁, s₂]
      end
    end
    nothing
  end
end

function _rmevent!(𝒬, s₁, s₂)
  intersection(Segment(s₁), Segment(s₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      BinaryTrees.delete!(𝒬, get(I))
    end
    nothing
  end
end

_isempty(𝒬) = isnothing(BinaryTrees.root(𝒬))

function _segmentswith(𝒟, p)
  P = typeof(p)
  S = Tuple{P,P}
  get(𝒟, p, S[])
end
