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
sample(object, method) = sample(Random.GLOBAL_RNG, object, method)

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

function sample(rng::AbstractRNG, geometry::Union{Multi,Polygon},
                method::ContinuousSamplingMethod)
  mesh = discretize(geometry, FIST(rng))
  sample(rng, mesh, method)
end

function sample(rng::AbstractRNG, geometry::Ngon,
                method::ContinuousSamplingMethod)
  mesh = discretize(geometry, Dehn1899())
  sample(rng, mesh, method)
end

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
    sample([rng], object, nsamples, [weights], replace=false)

Generate `nsamples` samples from spatial `object`
uniformly or using `weights`, with or without
replacement depending on `replace` option.
"""
sample(object::DomainOrData, nsamples::Int,
       weights::AbstractVector=[]; replace=false) =
  sample(Random.GLOBAL_RNG, object, nsamples, weights, replace)

function sample(rng::AbstractRNG,
                object::DomainOrData,
                nsamples::Int,
                weights::AbstractVector,
                replace::Bool)
  method = if isempty(weights)
    UniformSampling(nsamples, replace)
  else
    WeightedSampling(nsamples, weights, replace)
  end
  sample(rng, object, method)
end
