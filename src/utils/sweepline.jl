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
  𝒬 = AVLTree{P}()
  𝒯 = AVLTree{S}()
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
  ss = ℒ[p] ∪ 𝒰[p] ∪ 𝒞[p]
  if length(ss) > 1
    I[p] = ss
  end
  for s in ℒ[p] ∪ 𝒞[p]
    BinaryTrees.delete!(𝒯, s)
  end
  for s in 𝒰[p] ∪ 𝒞[p]
    BinaryTrees.insert!(𝒯, s)
  end
  if length(𝒰[p] ∪ 𝒞[p]) == 0
    # n = BinaryTrees.search(𝒯, p)
  else
  end
end

_key(node::BinaryTrees.AVLNode) = node.key