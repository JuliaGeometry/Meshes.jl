# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ProductPartition(p₁, p₂)

A method for partitioning spatial objects using the product of two
partitioning methods `p₁` and `p₂`.
"""
struct ProductPartition{P1,P2} <: PartitionMethod
  p₁::P1
  p₂::P2
end

# general case
function partition(object, method::ProductPartition)
  # individual partition results
  s₁ = indices(partition(object, method.p₁))
  s₂ = indices(partition(object, method.p₂))

  # label-based representation
  l₁ = Vector{Int}(undef, nelements(object))
  l₂ = Vector{Int}(undef, nelements(object))
  for (i, inds) in enumerate(s₁)
    l₁[inds] .= i
  end
  for (i, inds) in enumerate(s₂)
    l₂[inds] .= i
  end

  # product of labels
  counter = 0
  d = Dict{Tuple{Int,Int},Int}()
  l = Vector{Int}(undef, nelements(object))
  for i in 1:nelements(object)
    t = (l₁[i], l₂[i])
    if t ∉ keys(d)
      counter += 1
      d[t] = counter
    end
    l[i] = d[t]
  end

  # return partition using label predicate
  pred(i, j) = l[i] == l[j]
  partition(object, PredicatePartition(pred))
end

# predicate partition method
function partition(object, method::ProductPartition{P1,P2}) where {P1<:PredicatePartitionMethod,
                                                                   P2<:PredicatePartitionMethod}
  pred(i, j) = method.p₁(i, j) * method.p₂(i, j)
  partition(object, PredicatePartition(pred))
end

# spatial predicate partition method
function partition(object, method::ProductPartition{P1,P2}) where {P1<:SPredicatePartitionMethod,
                                                                   P2<:SPredicatePartitionMethod}
  pred(x, y) = method.p₁(x, y) * method.p₂(x, y)
  partition(object, SpatialPredicatePartition(pred))
end