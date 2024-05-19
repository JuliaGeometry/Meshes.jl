# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallPartition(radius; metric=Euclidean())

A method for partitioning spatial objects into balls of a given
`radius` using a `metric`.
"""
struct BallPartition{ℒ<:Len,M} <: SPredicatePartitionMethod
  radius::ℒ
  metric::M
  BallPartition(radius::ℒ, metric::M) where {ℒ<:Len,M} = new{float(ℒ),M}(radius, metric)
end

BallPartition(radius, metric) = BallPartition(addunit(radius, u"m"), metric)

BallPartition(radius; metric=Euclidean()) = BallPartition(radius, metric)

(p::BallPartition)(x, y) = evaluate(p.metric, x, y) < p.radius
