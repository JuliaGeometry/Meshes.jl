# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformSampling(size, replace=false, ordered=false)

Sample elements uniformly from a given domain/data. Produce a
sample of given `size` with or without replacement depending on
the `replace` option. The option `ordered` can be used to return
samples in the same order of the domain/data.
"""
struct UniformSampling <: DiscreteSamplingMethod
  size::Int
  replace::Bool
  ordered::Bool
end

UniformSampling(size; replace=false, ordered=false) = UniformSampling(size, replace, ordered)

function sampleinds(rng::AbstractRNG, d::Domain, method::UniformSampling)
  s = method.size
  r = method.replace
  o = method.ordered
  m = WeightedSampling(s; replace=r, ordered=o)
  sampleinds(rng, d, m)
end
