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
struct BallSampling{T,M} <: DiscreteSamplingMethod
  radius::T
  metric::M
  maxsize::Union{Int,Nothing}
end

BallSampling(radius; metric=Euclidean(), maxsize=nothing) =
  BallSampling(radius, metric, maxsize)

function sample(rng::AbstractRNG, object, method::BallSampling)
  radius = method.radius
  metric = method.metric
  msize  = method.maxsize ≠ nothing ? method.maxsize : Inf

  # neighborhood search with ball
  ball = IsotropicBall(radius, metric)
  searcher = BallSearch(object, ball)

  locations = Vector{Int}()
  notviewed = trues(nelements(object))
  while length(locations) < msize && any(notviewed)
    location = rand(rng, findall(notviewed))
    pₒ = centroid(object, location)

    # neighbors (including the location)
    neighbors = search(pₒ, searcher)

    push!(locations, location)
    notviewed[neighbors] .= false
  end

  view(object, locations)
end
