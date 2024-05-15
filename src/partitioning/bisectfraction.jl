# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectFractionPartition(normal, fraction=0.5, maxiter=10)

A method for partitioning spatial objects into two half spaces
defined by a `normal` direction and a `fraction` of points.
The partition is returned within `maxiter` bisection iterations.
"""
struct BisectFractionPartition{V<:Vec} <: PartitionMethod
  normal::V
  fraction::Float64
  maxiter::Int
  BisectFractionPartition{V}(normal, fraction, maxiter) where {V<:Vec} = new(unormalize(normal), fraction, maxiter)
end

BisectFractionPartition(normal::V, fraction=0.5, maxiter=10) where {V<:Vec} =
  BisectFractionPartition{V}(normal, fraction, maxiter)

BisectFractionPartition(normal::Tuple, fraction=0.5, maxiter=10) =
  BisectFractionPartition(Vec(normal), fraction, maxiter)

function partitioninds(rng::AbstractRNG, domain::Domain, method::BisectFractionPartition)
  u = unit(lentype(domain))
  bbox = boundingbox(domain)
  n = method.normal
  f = method.fraction
  c = coordinates(center(bbox))
  d = diagonal(bbox)

  # maximum number of bisections
  maxiter = method.maxiter

  iter = 0
  a = c - d / 2u * n
  b = c + d / 2u * n
  subsets = Vector{Int}[]
  metadata = Dict()
  while iter < maxiter
    m = (a + b) / 2

    bisectpoint = BisectPointPartition(n, Point(m))
    subsets, metadata = partitioninds(rng, domain, bisectpoint)

    g = length(subsets[1]) / nelements(domain)

    g ≈ f && break
    g > f && (b = m)
    g < f && (a = m)

    iter += 1
  end

  subsets, metadata
end
