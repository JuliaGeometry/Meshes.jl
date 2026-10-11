# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallPartition(radius; distance=EuclideanDistance())

A method for partitioning objects into balls of a given `radius` using a `distance`.
"""
struct BallPartition{ℒ<:Len,D<:GeometricDistance} <: PointPredicatePartitionMethod
  radius::ℒ
  distance::D
  BallPartition(radius::ℒ, distance::D) where {ℒ<:Len,D} = new{float(ℒ),D}(radius, distance)
end

BallPartition(radius, distance) = BallPartition(aslen(radius), distance)

BallPartition(radius; distance=EuclideanDistance()) = BallPartition(radius, distance)

(p::BallPartition)(pᵢ, pⱼ) = p.distance(pᵢ, pⱼ) < p.radius
