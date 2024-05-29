# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformPartition(k, [shuffle])

A method for partitioning spatial objects uniformly into `k` subsets
of approximately equal size. Optionally `shuffle` the data (default
to `true`).
"""
struct UniformPartition <: PartitionMethod
  k::Int
  shuffle::Bool
end

UniformPartition(k) = UniformPartition(k, true)

function partitioninds(rng::AbstractRNG, domain::Domain, method::UniformPartition)
  n = nelements(domain)
  k = method.k

  assertion(k ≤ n, "number of subsets must be smaller than number of points")

  inds = method.shuffle ? shuffle(rng, 1:n) : collect(1:n)
  subsets = collect(Iterators.partition(inds, n ÷ k))

  subsets, Dict()
end
