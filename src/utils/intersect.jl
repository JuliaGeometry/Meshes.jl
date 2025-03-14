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
  # adjust vertices of segments
  segs = map(segments) do s
    a, b = extrema(s)
    a > b ? reverse(s) : s
  end

  # retrieve types
  s = first(segs)
  p = minimum(s)
  P = typeof(p)
  S = Tuple{P,P}

  # initialization
  𝒬 = BinaryTrees.AVLTree{P}()
  ℛ = BinaryTrees.AVLTree{S}()
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
    lookup[(a, b)] = i
  end

  # sweep line
  points = Vector{P}()
  seginds = Vector{Vector{Int}}()
  while !isnothing(BinaryTrees.root(𝒬))
    p = BinaryTrees.key(BinaryTrees.minnode(𝒬))
    BinaryTrees.delete!(𝒬, p)
    _handle!(points, seginds, lookup, p, S, 𝒬, ℛ, ℒ, 𝒰, 𝒞)
  end
  points, seginds
end

function _handle!(points, seginds, lookup, p, S, 𝒬, ℛ, ℒ, 𝒰, 𝒞)
  ℬ = get(ℒ, p, S[])
  ℰ = get(𝒰, p, S[])
  ℐ = get(𝒞, p, S[])
  _processend!(ℰ, 𝒬, ℛ, 𝒞)
  _processbegin!(ℬ, 𝒬, ℛ, 𝒞)
  _processintersects!(ℐ, 𝒬, ℛ, 𝒞)
  segs = ℬ ∪ ℰ ∪ ℐ
  inds = [lookup[s] for s in segs]
  if !isempty(segs)
    push!(points, p)
    push!(seginds, inds)
  end
end

function _processbegin!(ℬ, 𝒬, ℛ, 𝒞)
  for s in ℬ
    BinaryTrees.insert!(ℛ, s)
  end
  for s in ℬ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev) && !isnothing(next)
      _newevent!(𝒬, 𝒞, BinaryTrees.key(next), BinaryTrees.key(prev))
    end
    if !isnothing(prev)
      _newevent!(𝒬, 𝒞, BinaryTrees.key(prev), s)
    end
    if !isnothing(next)
      _newevent!(𝒬, 𝒞, s, BinaryTrees.key(next))
    end
  end
end

function _processend!(ℰ, 𝒬, ℛ, 𝒞)
  for s in ℰ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev) && !isnothing(next)
      _newevent!(𝒬, 𝒞, BinaryTrees.key(next), BinaryTrees.key(prev))
    end
    BinaryTrees.delete!(ℛ, s)
  end
end

function _processintersects!(ℐ, 𝒬, ℛ, 𝒞)
  for s in ℐ
    prev, _ = BinaryTrees.prevnext(ℛ, s)
    if !isnothing(prev)
      # find segments r and u
      r, _ = BinaryTrees.prevnext(ℛ, s)
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
        _newevent!(𝒬, 𝒞, BinaryTrees.key(r), BinaryTrees.key(prev))
      end
      if !isnothing(u)
        _newevent!(𝒬, 𝒞, BinaryTrees.key(u), s)
      end
    end
  end
end

function _newevent!(𝒬, 𝒞, (a₁, b₁), (a₂, b₂))
  intersection(Segment(a₁, b₁), Segment(a₂, b₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      p = get(I)
      BinaryTrees.insert!(𝒬, p)
      if haskey(𝒞, p)
        push!(𝒞[p], (a₁, b₁), (a₂, b₂))
      else
        𝒞[p] = [(a₁, b₁), (a₂, b₂)]
      end
    end
    nothing
  end
end

function _rmevent!(𝒬, (a₁, b₁), (a₂, b₂))
  intersection(Segment(a₁, b₁), Segment(a₂, b₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      BinaryTrees.delete!(𝒬, get(I))
    end
    nothing
  end
end
