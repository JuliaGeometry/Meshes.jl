# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ProductPartition(first, second)

A method for partitioning objects using the product of
`first` and `second` partitioning methods.
"""
struct ProductPartition{P1,P2} <: PartitionMethod
  first::P1
  second::P2
end

# general case
function partitioninds(rng::AbstractRNG, domain::Domain, method::ProductPartition)
  # individual partition results
  s₁, _ = partitioninds(rng, domain, method.first)
  s₂, _ = partitioninds(rng, domain, method.second)

  # label-based representation
  l₁ = Vector{Int}(undef, nelements(domain))
  l₂ = Vector{Int}(undef, nelements(domain))
  for (i, inds) in enumerate(s₁)
    l₁[inds] .= i
  end
  for (i, inds) in enumerate(s₂)
    l₂[inds] .= i
  end

  # product of labels
  counter = 0
  d = Dict{Tuple{Int,Int},Int}()
  l = Vector{Int}(undef, nelements(domain))
  for i in 1:nelements(domain)
    t = (l₁[i], l₂[i])
    if t ∉ keys(d)
      counter += 1
      d[t] = counter
    end
    l[i] = d[t]
  end

  # return partition using label predicate
  pred(i, j) = l[i] == l[j]
  partitioninds(rng, domain, IndexPredicatePartition(pred))
end

# index predicate partition method
function partitioninds(
  rng::AbstractRNG,
  domain::Domain,
  method::ProductPartition{P1,P2}
) where {P1<:IndexPredicatePartitionMethod,P2<:IndexPredicatePartitionMethod}
  pred(i, j) = method.first(i, j) * method.second(i, j)
  partitioninds(rng, domain, IndexPredicatePartition(pred))
end

# point predicate partition method
function partitioninds(
  rng::AbstractRNG,
  domain::Domain,
  method::ProductPartition{P1,P2}
) where {P1<:PointPredicatePartitionMethod,P2<:PointPredicatePartitionMethod}
  pred(pᵢ, pⱼ) = method.first(pᵢ, pⱼ) * method.second(pᵢ, pⱼ)
  partitioninds(rng, domain, PointPredicatePartition(pred))
end
