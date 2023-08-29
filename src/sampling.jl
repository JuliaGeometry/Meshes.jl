# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SamplingMethod

A method for sampling from geometric objects.
"""
abstract type SamplingMethod end

"""
    sample([rng], object, method)

Sample elements or points from geometric `object`
with `method`. Optionally, specify random number
generator `rng`.
"""
sample(object, method::SamplingMethod) = sample(Random.GLOBAL_RNG, object, method)

"""
    sampleinds(rng, domain, method)

Sample indices of elements in `domain` with discrete
sampling `method` and random number generator `rng`.
"""
function sampleinds end

"""
    DiscreteSamplingMethod

A method for sampling from discrete representations
of geometric objects such as meshes or collections
of geometries.
"""
abstract type DiscreteSamplingMethod <: SamplingMethod end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sampling/uniform.jl")
include("sampling/weighted.jl")
include("sampling/ball.jl")
include("sampling/block.jl")

"""
    ContinuousSamplingMethod

A method for sampling from continuous representations
of geometric objects. In this case, geometric objects
are interpreted as a set of points in the embedding
space.
"""
abstract type ContinuousSamplingMethod <: SamplingMethod end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sampling/regular.jl")
include("sampling/homogeneous.jl")
include("sampling/mindistance.jl")

# ----------
# FALLBACKS
# ----------

sample(rng::AbstractRNG, d::Domain, method::DiscreteSamplingMethod) = view(d, sampleinds(rng, d, method))

sample(rng::AbstractRNG, g::Geometry, method::ContinuousSamplingMethod) = sample(rng, discretize(g), method)

# ----------
# UTILITIES
# ----------

"""
    sample([rng], domain, size, [weights]; replace=false, ordered=false)

Generate `size` samples from `domain` uniformly or using `weights`,
with or without replacement depending on the `replace` option. The
option `ordered` can be used to return samples in the same order of
the `domain`.
"""
function sample(domain::Domain, size::Int, weights=nothing; replace=false, ordered=false)
  method = WeightedSampling(size, weights; replace=replace, ordered=ordered)
  sample(Random.GLOBAL_RNG, domain, method)
end
