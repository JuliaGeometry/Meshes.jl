# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DirectionPartition(direction; [tol])

A method for partitioning spatial objects along a given `direction`
with bandwidth tolerance `tol`.
"""
struct DirectionPartition{V<:Vec,ℒ<:Len} <: SPredicatePartitionMethod
  direction::V
  tol::ℒ
  DirectionPartition(direction::V, tol::ℒ) where {V<:Vec,ℒ<:Len} = new{V,float(ℒ)}(unormalize(direction), tol)
end

DirectionPartition(direction::Vec, tol) = DirectionPartition(direction, addunit(tol, u"m"))

DirectionPartition(direction::Vec; tol=atol(eltype(direction))) = DirectionPartition(direction, tol)

DirectionPartition(direction::Tuple; kwargs...) = DirectionPartition(Vec(direction); kwargs...)

function (p::DirectionPartition)(x, y)
  δ = x - y
  d = p.direction
  k = ustrip(δ ⋅ d)
  norm(δ - k * d) < p.tol
end
