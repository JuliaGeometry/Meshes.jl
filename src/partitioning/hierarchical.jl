# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HierarchicalPartition(first, second)

A partitioning method in which a `first` partition is applied
and then a `second` partition is applied to each subset of the
`first`.
"""
struct HierarchicalPartition <: PartitionMethod
  first::PartitionMethod
  second::PartitionMethod
end

function partitioninds(rng::AbstractRNG, domain::Domain, method::HierarchicalPartition)
  subsets = Vector{Int}[]

  # use first partition method
  s₁, _ = partitioninds(rng, domain, method.first)

  # use second method to partition the first
  for is in s₁
    v = view(domain, is)
    s₂, _ = partitioninds(rng, v, method.second)
    for js in s₂
      push!(subsets, is[js])
    end
  end

  subsets, Dict()
end
