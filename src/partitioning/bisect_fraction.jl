# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectFractionPartition(normal, fraction=0.5, maxiter=10)

A method for partitioning spatial data into two half spaces
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

function partition(object, method::BisectFractionPartition)
  bbox = boundingbox(object)
  n = method.normal
  f = method.fraction
  c = coordinates(center(bbox))
  d = diagonal(bbox)

  # maximum number of bisections
  maxiter = method.maxiter

  iter = 0; p = 0
  a = c - d/2 * n
  b = c + d/2 * n
  while iter < maxiter
    m = (a + b) / 2

    p = partition(object, BisectPointPartition(n, Point(m)))
    g = nelements(p[1]) / nelements(object)

    g ≈ f && break
    g > f && (b = m)
    g < f && (a = m)

    iter += 1
  end

  p
end
