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

sample(rng::AbstractRNG, geom::Geometry, method::ContinuousSamplingMethod) =
  sample(rng, discretize(geom), method)

# ----------------
# IMPLEMENTATIONS
# ----------------

# discrete sampling
include("sampling/uniform.jl")
include("sampling/weighted.jl")
include("sampling/ball.jl")
include("sampling/block.jl")

# continuous sampling
include("sampling/regular.jl")
include("sampling/homogeneous.jl")
include("sampling/mindistance.jl")

# ----------
# UTILITIES
# ----------

"""
    sample([rng], object, size, [weights]; replace=false, ordered=false)

Generate `size` samples from `object` uniformly or using `weights`,
with or without replacement depending on the `replace` option. The
option `ordered` can be used to return samples in the same order of
the `object`.
"""
function sample(
  object::DomainOrData,
  size::Int,
  weights=nothing;
  replace=false,
  ordered=false
)
  method = WeightedSampling(size, weights; replace=replace, ordered=ordered)
  sample(Random.GLOBAL_RNG, object, method)
end
