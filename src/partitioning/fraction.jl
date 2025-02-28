# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FractionPartition(fraction; shuffle=true)

A method for partitioning objects according to a given `fraction`.
Optionally `shuffle` elements before partitioning.
"""
struct FractionPartition <: PartitionMethod
  fraction::Float64
  shuffle::Bool
end

function FractionPartition(fraction; shuffle=true)
  assertion(0 < fraction < 1, "fraction must be in interval (0,1)")
  FractionPartition(fraction, shuffle)
end

function partitioninds(rng::AbstractRNG, domain::Domain, method::FractionPartition)
  n = nelements(domain)
  f = round(Int, method.fraction * n)

  locs = method.shuffle ? randperm(rng, n) : 1:n
  subsets = [locs[1:f], locs[(f + 1):n]]

  subsets, Dict()
end
