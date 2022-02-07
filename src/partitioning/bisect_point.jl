# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectPointPartition(normal, point)

A method for partitioning spatial objects into two half spaces
defined by a `normal` direction and a reference `point`.
"""
struct BisectPointPartition{Dim,T} <: PartitionMethod
  normal::Vec{Dim,T}
  point::Point{Dim,T}

  function BisectPointPartition{Dim,T}(normal, point) where {Dim,T}
    new(normalize(normal), point)
  end
end

BisectPointPartition(normal::Vec{Dim,T}, point::Point{Dim,T}) where {Dim,T} =
  BisectPointPartition{Dim,T}(normal, point)

BisectPointPartition(normal::NTuple{Dim,T}, point::NTuple{Dim,T}) where {Dim,T} =
 BisectPointPartition(Vec(normal), Point(point))

function partition(::AbstractRNG, object, method::BisectPointPartition)
  n = method.normal
  p = method.point

  left, right = Int[], Int[]
  for location in 1:nelements(object)
    pₒ = centroid(object, location)
    if (pₒ - p) ⋅ n < zero(coordtype(object))
      push!(left, location)
    else
      push!(right, location)
    end
  end

  Partition(object, [left,right])
end
