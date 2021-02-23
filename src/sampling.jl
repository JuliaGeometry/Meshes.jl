# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SamplingMethod

A method for sampling from geometries.
"""
abstract type SamplingMethod end

"""
    sample(geometry, method)

Sample elements from `geometry` with `method`.
"""
function sample end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sampling/regular.jl")
include("sampling/uniform.jl")
include("sampling/weighted.jl")
include("sampling/ball.jl")

# ----------
# UTILITIES
# ----------

"""
    sample(object, nsamples, [weights], replace=false)

Generate `nsamples` samples from spatial `object`
uniformly or using `weights`, with or without
replacement depending on `replace` option.
"""
function sample(object::Union{Domain,Data}, nsamples::Int,
                weights::AbstractVector=[]; replace=false)
  if isempty(weights)
    sample(object, UniformSampling(nsamples, replace))
  else
    sample(object, WeightedSampling(nsamples, weights, replace))
  end
end
