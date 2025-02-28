# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DirectionPartition(direction; [tol])

A method for partitioning objects along a given `direction` with bandwidth tolerance `tol`.
"""
struct DirectionPartition{V<:Vec,ℒ<:Len} <: PointPredicatePartitionMethod
  dir::V
  tol::ℒ
  DirectionPartition(dir::V, tol::ℒ) where {V<:Vec,ℒ<:Len} = new{V,float(ℒ)}(unormalize(dir), tol)
end

DirectionPartition(dir::Vec, tol) = DirectionPartition(dir, addunit(tol, u"m"))

DirectionPartition(dir::Vec; tol=atol(eltype(dir))) = DirectionPartition(dir, tol)

DirectionPartition(dir::Tuple; kwargs...) = DirectionPartition(Vec(dir); kwargs...)

function (p::DirectionPartition)(pᵢ, pⱼ)
  δ = pᵢ - pⱼ
  d = p.dir
  k = ustrip(δ ⋅ d)
  norm(δ - k * d) < p.tol
end
