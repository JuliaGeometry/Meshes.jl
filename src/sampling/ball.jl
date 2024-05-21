# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSampling(radius; [options])

A method for sampling isolated elements from a given domain/data
according to a norm-ball of given `radius`.

## Options

* `metric`  - Metric for the ball (default to `Euclidean()`)
* `maxsize` - Maximum size of the resulting sample (default to none)
"""
struct BallSampling{ℒ<:Len,M} <: DiscreteSamplingMethod
  radius::ℒ
  metric::M
  maxsize::Union{Int,Nothing}
  BallSampling(radius::ℒ, metric::M, maxsize) where {ℒ<:Len,M} = new{float(ℒ),M}(radius, metric, maxsize)
end

BallSampling(radius::Len; metric=Euclidean(), maxsize=nothing) = BallSampling(radius, metric, maxsize)

BallSampling(radius; kwargs...) = BallSampling(addunit(radius, u"m"); kwargs...)

function sampleinds(rng::AbstractRNG, d::Domain, method::BallSampling)
  radius = method.radius
  metric = method.metric
  msize = isnothing(method.maxsize) ? Inf : method.maxsize

  # neighborhood search with ball
  ball = MetricBall(radius, metric)
  searcher = BallSearch(d, ball)

  inds = Int[]
  notviewed = trues(nelements(d))
  while length(inds) < msize && any(notviewed)
    ind = rand(rng, findall(notviewed))
    pₒ = centroid(d, ind)

    # neighbors (including the index)
    neighbors = search(pₒ, searcher)

    push!(inds, ind)
    notviewed[neighbors] .= false
  end

  inds
end
