# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformPartition(k; shuffle=true)

A method for partitioning objects uniformly into `k` subsets
of approximately equal size. Optionally `shuffle` the data.
"""
struct UniformPartition <: PartitionMethod
  k::Int
  shuffle::Bool
end

UniformPartition(k; shuffle=true) = UniformPartition(k, shuffle)

function partitioninds(rng::AbstractRNG, domain::Domain, method::UniformPartition)
  n = nelements(domain)
  k = method.k

  assertion(k ≤ n, "number of subsets must be smaller than number of points")

  inds = method.shuffle ? shuffle(rng, 1:n) : collect(1:n)
  subsets = collect(Iterators.partition(inds, n ÷ k))

  subsets, Dict()
end
