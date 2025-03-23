# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    bentleyottmann(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O(n⋅log(n)) time using
Bentley-Ottmann sweep line algorithm.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function bentleyottmann(segments; digits=_digits(segments))
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
  counter = 1
  points = Vector{P}()
  seginds = Vector{Vector{Int}}()
  visited = Dict{P,Int}()
  while !BinaryTrees.isempty(𝒬)
    # current point (or event)
    p = BinaryTrees.key(BinaryTrees.minnode(𝒬))

    # delete point from event queue
    BinaryTrees.delete!(𝒬, p)

    # handle event, i.e. update 𝒬, ℛ and ℳ
    ℬₚ = get(ℬ, p, S[]) # segments with p at the begin
    ℰₚ = get(ℰ, p, S[]) # segments with p at the end
    ℳₚ = get(ℳ, p, S[]) # segments with p at the middle
    _handlebeg!(ℬₚ, 𝒬, ℛ, ℳ, digits)
    _handleend!(ℰₚ, 𝒬, ℛ, ℳ, digits)
    _handlemid!(ℳₚ, 𝒬, ℛ, ℳ, digits)

    # report intersection point and segment indices
    inds = [lookup[s] for s in ℳₚ]
    if !isempty(inds)
      if p ∈ keys(visited)
        seginds[visited[p]] = inds
      else
        push!(points, p)
        push!(seginds, inds)
        push!(visited, p => counter)
        counter += 1
      end
    end
  end

  points, seginds
end

function _handlebeg!(ℬₚ, 𝒬, ℛ, ℳ, digits)
  for s in ℬₚ
    BinaryTrees.insert!(ℛ, s)
    prev, next = BinaryTrees.prevnext(ℛ, s)
    isnothing(prev) || _newevent!(𝒬, ℳ, BinaryTrees.key(prev), s, digits)
    isnothing(next) || _newevent!(𝒬, ℳ, s, BinaryTrees.key(next), digits)
    isnothing(prev) || isnothing(next) || _rmevent!(𝒬, s, s, digits)
  end
end

function _handleend!(ℰₚ, 𝒬, ℛ, ℳ, digits)
  for s in ℰₚ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    isnothing(prev) || isnothing(next) || _newevent!(𝒬, ℳ, BinaryTrees.key(prev), BinaryTrees.key(next), digits)
    BinaryTrees.delete!(ℛ, s)
  end
end

function _handlemid!(ℳₚ, 𝒬, ℛ, ℳ, digits)
  for s in ℳₚ
    prev, next = BinaryTrees.prevnext(ℛ, s)
    r = !isnothing(prev) ? BinaryTrees.key(prev) : nothing
    t = !isnothing(next) ? BinaryTrees.key(next) : nothing
    if !isnothing(r)
      _newevent!(𝒬, ℳ, r, s, digits)
      if !isnothing(t)
        _newevent!(𝒬, ℳ, r, t, digits)
      end
    end
    if !isnothing(t)
      _, next = BinaryTrees.prevnext(ℛ, BinaryTrees.key(next))
      u = !isnothing(next) ? BinaryTrees.key(next) : nothing
      if !isnothing(u)
        _newevent!(𝒬, ℳ, t, u, digits)
        if !isnothing(r)
          _newevent!(𝒬, ℳ, r, u, digits)
        end
      end
    end
  end
end

function _newevent!(𝒬, ℳ, s₁, s₂, digits)
  intersection(Segment(s₁), Segment(s₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      p = roundcoords(get(I); digits)
      if haskey(ℳ, p)
        if s₁ ∉ ℳ[p]
          push!(ℳ[p], s₁)
          BinaryTrees.insert!(𝒬, p)
        end
        if s₂ ∉ ℳ[p]
          push!(ℳ[p], s₂)
        end
      else
        ℳ[p] = [s₁, s₂]
        BinaryTrees.insert!(𝒬, p)
      end
    end
    nothing
  end
end

function _rmevent!(𝒬, s₁, s₂, digits)
  intersection(Segment(s₁), Segment(s₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      p = roundcoords(get(I); digits)
      BinaryTrees.delete!(𝒬, p)
    end
    nothing
  end
end

function _digits(segments)
  s = first(segments)
  ℒ = lentype(s)
  τ = ustrip(atol(ℒ))
  round(Int, -log10(τ)) - 1
end
