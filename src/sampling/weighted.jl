# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WeightedSampling(size, [weights]; replace=false)

Sample elements from a given domain/data using `weights`. Produce a sample
of given `size` with or without replacement depending on the `replace`
option. By default weights are uniform.
"""
struct WeightedSampling{W} <: DiscreteSamplingMethod
  size::Int
  weights::W
  replace::Bool
end

WeightedSampling(size, weights=nothing; replace=false) =
  WeightedSampling(size, weights, replace)

function sample(rng::AbstractRNG, Ω::DomainOrData, method::WeightedSampling)
  n = nelements(Ω)
  s = method.size
  w = method.weights
  r = method.replace
  if s > n && r == false
    @error "invalid sample size for sampling without replacement"
  end
  ws = isnothing(w) ? fill(1/n, n) : collect(w)
  @assert length(ws) == n "invalid number of weights for object"
  view(Ω, sample(rng, 1:n, Weights(ws), s, replace=r))
end
