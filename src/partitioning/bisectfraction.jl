# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectFractionPartition(normal, fraction=0.5, maxiter=10)

A method for partitioning spatial objects into two half spaces
defined by a `normal` direction and a `fraction` of points.
The partition is returned within `maxiter` bisection iterations.
"""
struct BisectFractionPartition{Dim,T} <: PartitionMethod
  normal::Vec{Dim,T}
  fraction::Float64
  maxiter::Int

  function BisectFractionPartition{Dim,T}(normal, fraction, maxiter) where {Dim,T}
    new(normalize(normal), fraction, maxiter)
  end
end

BisectFractionPartition(normal::Vec{Dim,T}, fraction=0.5, maxiter=10) where {Dim,T} =
  BisectFractionPartition{Dim,T}(normal, fraction, maxiter)

BisectFractionPartition(normal::NTuple{Dim,T}, fraction=0.5, maxiter=10) where {Dim,T} =
  BisectFractionPartition(Vec(normal), fraction, maxiter)

function partsubsets(rng::AbstractRNG, domain::Domain, method::BisectFractionPartition)
  bbox = boundingbox(domain)
  n = method.normal
  f = method.fraction
  c = coordinates(center(bbox))
  d = diagonal(bbox)

  # maximum number of bisections
  maxiter = method.maxiter

  iter = 0
  a = c - d / 2 * n
  b = c + d / 2 * n
  subsets = Vector{Int}[]
  metadata = Dict()
  while iter < maxiter
    m = (a + b) / 2

    bisectpoint = BisectPointPartition(n, Point(m))
    subsets, metadata = partsubsets(rng, domain, bisectpoint)

    g = length(subsets[1]) / nelements(domain)

    g ≈ f && break
    g > f && (b = m)
    g < f && (a = m)

    iter += 1
  end

  subsets, metadata
end
