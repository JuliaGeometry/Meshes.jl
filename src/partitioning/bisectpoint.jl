# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BisectPointPartition(normal, point)

A method for partitioning spatial objects into two half spaces
defined by a `normal` direction and a reference `point`.
"""
struct BisectPointPartition{Dim,V<:Vec{Dim},P<:Point{Dim}} <: PartitionMethod
  normal::V
  point::P
  BisectPointPartition{Dim,V,P}(normal, point) where {Dim,V<:Vec{Dim},P<:Point{Dim}} = new(unormalize(normal), point)
end

BisectPointPartition(normal::V, point::P) where {Dim,V<:Vec{Dim},P<:Point{Dim}} =
  BisectPointPartition{Dim,V,P}(normal, point)

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
