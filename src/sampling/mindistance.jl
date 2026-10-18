# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MinDistanceSampling(α, ρ=0.65, δ=100, metric=EuclideanDistance())

Generate sample from geometric object such that all pairs of
points are at least `α` units of distance away from each other.
Optionally specify the relative radius `ρ` for the packing
pattern, the oversampling factor `δ` and the `metric`.

This method is sometimes referred to as Poisson disk sampling
or blue noise sampling in the computer graphics community.

## References

* Lagae, A. & Dutré, P. 2007. [A Comparison of Methods for
  Generating Poisson Disk Distributions]
  (https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-8659.2007.01100.x)
* Bowers et al. 2010. [Parallel Poisson disk sampling with
  spectrum analysis on surfaces](https://dl.acm.org/doi/10.1145/1882261.1866188)
* Medeiros et al. 2014. [Fast adaptive blue noise on polygonal surfaces]
  (https://www.sciencedirect.com/science/article/abs/pii/S1524070313000313)
"""
struct MinDistanceSampling{ℒ<:Len,M} <: ContinuousSamplingMethod
  α::ℒ
  ρ::ℒ
  δ::Int
  metric::M
  MinDistanceSampling(α::ℒ, ρ::ℒ, δ, metric::M) where {ℒ<:Len,M} = new{float(ℒ),M}(α, ρ, δ, metric)
end

MinDistanceSampling(α::Len, ρ::Len, δ, metric) = MinDistanceSampling(promote(α, ρ)..., δ, metric)

MinDistanceSampling(α, ρ, δ, metric) = MinDistanceSampling(aslen(α), aslen(ρ), δ, metric)

MinDistanceSampling(α::T; ρ=T(0.65), δ=100, metric=EuclideanDistance()) where {T} = MinDistanceSampling(α, ρ, δ, metric)

sample(rng::AbstractRNG, d::Domain, method::MinDistanceSampling) = _sample(rng, d, method)

sample(rng::AbstractRNG, b::Ball, method::MinDistanceSampling) = _sample(rng, b, method)

function _sample(rng::AbstractRNG, obj, method::MinDistanceSampling)
  # retrieve parameters
  α = method.α
  ρ = method.ρ
  δ = method.δ
  m = method.metric

  # total volume/area of the object
  V = measure(obj)

  # expected number of Poisson samples
  # for relative radius (Lagae & Dutré 2007)
  N = 2V / √3 * (ρ / α)^2

  # number of oversamples (Medeiros et al. 2014)
  O = ceil(Int, δ * ustrip(N))

  # oversample the object
  points = sample(rng, obj, HomogeneousSampling(O))

  # collect points into point set
  𝒫 = PointSet(collect(points))

  # discard points that do not satisfy distance criterion
  sample(rng, 𝒫, BallSampling(α, metric=m))
end
