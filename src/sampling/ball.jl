# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSampling(radius; distance=EuclideanDistance(), maxsize=nothing)

A method for sampling isolated elements from a given domain/data
according to a norm-ball of given `radius`. Optionally, specify
the `distance` used for the metric ball and the `maxsize` of the
resulting sample.
"""
struct BallSampling{ℒ<:Len,D} <: DiscreteSamplingMethod
  radius::ℒ
  distance::D
  maxsize::Union{Int,Nothing}
  BallSampling(radius::ℒ, distance::D, maxsize) where {ℒ<:Len,D} = new{float(ℒ),D}(radius, distance, maxsize)
end

BallSampling(radius::Len; distance=EuclideanDistance(), maxsize=nothing) = BallSampling(radius, distance, maxsize)

BallSampling(radius; kwargs...) = BallSampling(aslen(radius); kwargs...)

function sampleinds(rng::AbstractRNG, d::Domain, method::BallSampling)
  radius = method.radius
  distance = method.distance
  msize = isnothing(method.maxsize) ? Inf : method.maxsize

  # neighborhood search with ball
  ball = MetricBall(radius, distance)
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
