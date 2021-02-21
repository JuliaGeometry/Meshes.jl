# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformSampling(size, replace=false)

Sample elements uniformly from a given domain. Produce a sample
of given `size` with or without replacement depending on the
`replace` option.
"""
struct UniformSampling <: SamplingMethod
  size::Int
  replace::Bool
end

UniformSampling(size::Int) = UniformSampling(size, false)

function sample(domain::Domain, method::UniformSampling)
  n = nelements(domain)
  s = method.size
  r = method.replace
  if s > n && r == false
    @error "invalid sample size for sampling without replacement"
  end
  view(domain, sample(1:n, s, replace=r))
end
