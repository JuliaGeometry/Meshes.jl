# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectPointPartition(normal, point)

A method for partitioning objects into two half spaces
defined by a `normal` direction and a reference `point`.
"""
struct BisectPointPartition{V<:Vec,P<:Point} <: PartitionMethod
  normal::V
  point::P
  BisectPointPartition{V,P}(normal, point) where {V<:Vec,P<:Point} = new(unormalize(normal), point)
end

BisectPointPartition(normal::V, point::P) where {V<:Vec,P<:Point} = BisectPointPartition{V,P}(normal, point)

BisectPointPartition(normal::NTuple{Dim}, point::NTuple{Dim}) where {Dim} =
  BisectPointPartition(Vec(normal), Point(point))

function partitioninds(::AbstractRNG, domain::Domain, method::BisectPointPartition)
  n = method.normal
  p = method.point

  left, right = Int[], Int[]
  for location in 1:nelements(domain)
    pₒ = centroid(domain, location)
    if isnegative((pₒ - p) ⋅ n)
      push!(left, location)
    else
      push!(right, location)
    end
  end

  [left, right], Dict()
end
