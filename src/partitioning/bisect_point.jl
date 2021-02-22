# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectPointPartition(normal, point)

A method for partitioning spatial data into two half spaces
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

function partition(object, method::BisectPointPartition)
  Dim = embeddim(object)
  T = coordtype(object)
  
  n = method.normal
  p = method.point

  x = MVector{Dim,T}(undef)

  left  = Vector{Int}()
  right = Vector{Int}()
  for location in 1:nelements(object)
    coordinates!(x, object, location)
    if (Point(x) - p) ⋅ n < zero(T)
      push!(left, location)
    else
      push!(right, location)
    end
  end

  Partition(object, [left,right])
end
