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
sample(object, method::SamplingMethod) = sample(Random.default_rng(), object, method)

"""
    DiscreteSamplingMethod

A method for sampling from discrete representations
of geometric objects such as meshes or collections
of geometries.
"""
abstract type DiscreteSamplingMethod <: SamplingMethod end

sample(rng::AbstractRNG, object, method::DiscreteSamplingMethod) = view(object, sampleinds(rng, object, method))

"""
    sampleinds(rng, domain, method)

Sample indices of elements in `domain` with discrete
sampling `method` and random number generator `rng`.
"""
function sampleinds end

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

sample(rng::AbstractRNG, g::Geometry, method::ContinuousSamplingMethod) = sample(rng, discretize(g), method)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sampling/regular.jl")
include("sampling/homogeneous.jl")
include("sampling/mindistance.jl")
include("sampling/fibonacci.jl")
include("sampling/adaptive.jl")

# ----------
# UTILITIES
# ----------

"""
    sample([rng], domain, size, [weights]; replace=false, ordered=false)

Utility method that calls the `sample` function using `WeightedSampling(size, weights; replace, ordered)`.
If `weights` is not defined, this is equivalent to using `UniformSampling(size; replace, ordered)`.
"""
sample(domain::Domain, size::Int, weights=nothing; kwargs...) =
  sample(Random.default_rng(), domain, size, weights; kwargs...)
sample(rng::AbstractRNG, domain::Domain, size::Int, weights=nothing; kwargs...) =
  sample(rng, domain, WeightedSampling(size, weights; kwargs...))
