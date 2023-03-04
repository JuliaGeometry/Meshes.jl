# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallPartition(radius; metric=Euclidean())

A method for partitioning spatial objects into balls of a given
`radius` using a `metric`.
"""
struct BallPartition{T,M} <: SPredicatePartitionMethod
  radius::T
  metric::M
end

function BallPartition(radius::T; metric::M=Euclidean()) where {T,M}
  return BallPartition{T,M}(radius, metric)
end

(p::BallPartition)(x, y) = evaluate(p.metric, x, y) < p.radius
