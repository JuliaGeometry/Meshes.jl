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

function partition(object, method::HierarchicalPartition)
  result = Vector{Int}[]

  # use first partition method
  p = partition(object, method.first)

  # use second method to partition the first
  s = indices(p)
  for (i, d) in Iterators.enumerate(p)
    q = partition(d, method.second)

    for js in indices(q)
      push!(result, s[i][js])
    end
  end

  Partition(object, result)
end