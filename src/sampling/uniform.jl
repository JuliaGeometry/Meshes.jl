# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformSampling(size, replace=false)

Sample elements uniformly from a given domain/data. Produce a
sample of given `size` with or without replacement depending on
the `replace` option.
"""
struct UniformSampling <: DiscreteSamplingMethod
  size::Int
  replace::Bool
end

UniformSampling(size::Int) = UniformSampling(size, false)

function sample(object, method::UniformSampling)
  n = nelements(object)
  s = method.size
  r = method.replace
  if s > n && r == false
    @error "invalid sample size for sampling without replacement"
  end
  view(object, sample(1:n, s, replace=r))
end
