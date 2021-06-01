# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SamplingMethod

A method for sampling from geometric objects.
"""
abstract type SamplingMethod end

"""
    DiscreteSamplingMethod

A method for sampling from discrete representations
of geometric objects such as meshes or collections
of geometries.
"""
abstract type DiscreteSamplingMethod end

"""
    ContinuousSamplingMethod

A method for sampling from continuous representations
of geometric objects. In this case, geometric objects
are interpreted as a set of points in the embedding
space.
"""
abstract type ContinuousSamplingMethod end

"""
    sample(object, method)

Sample elements from `object` with `method`.
"""
function sample end

# ----------------
# IMPLEMENTATIONS
# ----------------

# discrete sampling
include("sampling/uniform.jl")
include("sampling/weighted.jl")
include("sampling/ball.jl")

# continuous sampling
include("sampling/regular.jl")
include("sampling/homogeneous.jl")
include("sampling/mindistance.jl")

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
