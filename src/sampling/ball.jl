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
struct BallSampling{T,M} <: SamplingMethod
  radius::T
  metric::M
  maxsize::Union{Int,Nothing}
end

BallSampling(radius; metric=Euclidean(), maxsize=nothing) =
  BallSampling(radius, metric, maxsize)

function sample(object, method::BallSampling)
  Dim = embeddim(object)
  T = coordtype(object)
  radius = method.radius
  metric = method.metric
  msize  = method.maxsize ≠ nothing ? method.maxsize : Inf

  # neighborhood search with ball
  ball = NormBall(radius, metric)
  searcher = NeighborhoodSearch(object, ball)

  # pre-allocate memory for coordinates
  coords = MVector{Dim,T}(undef)

  locations = Vector{Int}()
  notviewed = trues(nelements(object))
  while length(locations) < msize && any(notviewed)
    location = rand(findall(notviewed))
    coordinates!(coords, object, location)

    # neighbors (including the location)
    neighbors = search(Point(coords), searcher)

    push!(locations, location)
    notviewed[neighbors] .= false
  end

  view(object, locations)
end
