# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WeightedSampling(size, [weights]; replace=false, ordered=false)

Sample elements from a given domain/data using `weights`. Produce a sample
of given `size` with or without replacement depending on the `replace`
option. The option `ordered` can be used to return samples in the same
order of the original domain/data. By default weights are uniform.
"""
struct WeightedSampling{W} <: DiscreteSamplingMethod
  size::Int
  weights::W
  replace::Bool
  ordered::Bool
end

function WeightedSampling(size; replace=false, ordered=false)
  return WeightedSampling(size, nothing, replace, ordered)
end

function WeightedSampling(size, weights::AbstractWeights; replace=false, ordered=false)
  return WeightedSampling(size, weights, replace, ordered)
end

function WeightedSampling(size, weights::Nothing; replace=false, ordered=false)
  return WeightedSampling(size, weights, replace, ordered)
end

function WeightedSampling(size, weights; replace=false, ordered=false)
  return WeightedSampling(size, Weights(collect(weights)), replace, ordered)
end

function sample(rng::AbstractRNG, Ω::DomainOrData, method::WeightedSampling)
  n = nitems(Ω)
  s = method.size
  w = method.weights
  r = method.replace
  o = method.ordered

  if s > n && r == false
    throw(ArgumentError("invalid sample size for sampling without replacement"))
  end

  inds = if isnothing(w)
    sample(rng, 1:n, s; replace=r, ordered=o)
  else
    sample(rng, 1:n, w, s; replace=r, ordered=o)
  end

  return view(Ω, inds)
end
