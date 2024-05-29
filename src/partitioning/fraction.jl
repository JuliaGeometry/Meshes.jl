# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FractionPartition(fraction, shuffle=true)

A method for partitioning spatial objects according to a given `fraction`.
Optionally `shuffle` elements before partitioning.
"""
struct FractionPartition <: PartitionMethod
  fraction::Float64
  shuffle::Bool

  function FractionPartition(fraction, shuffle)
    assertion(0 < fraction < 1, "fraction must be in interval (0,1)")
    new(fraction, shuffle)
  end
end

FractionPartition(fraction) = FractionPartition(fraction, true)

function partitioninds(rng::AbstractRNG, domain::Domain, method::FractionPartition)
  n = nelements(domain)
  f = round(Int, method.fraction * n)

  locs = method.shuffle ? randperm(rng, n) : 1:n
  subsets = [locs[1:f], locs[(f + 1):n]]

  subsets, Dict()
end
