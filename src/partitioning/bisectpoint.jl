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
    return new(normalize(normal), point)
  end
end

function BisectPointPartition(normal::Vec{Dim,T}, point::Point{Dim,T}) where {Dim,T}
  return BisectPointPartition{Dim,T}(normal, point)
end

function BisectPointPartition(normal::NTuple{Dim,T}, point::NTuple{Dim,T}) where {Dim,T}
  return BisectPointPartition(Vec(normal), Point(point))
end

function partsubsets(::AbstractRNG, domain::Domain, method::BisectPointPartition)
  n = method.normal
  p = method.point

  left, right = Int[], Int[]
  for location in 1:nelements(domain)
    pₒ = centroid(domain, location)
    if (pₒ - p) ⋅ n < zero(coordtype(domain))
      push!(left, location)
    else
      push!(right, location)
    end
  end

  return [left, right], Dict()
end
