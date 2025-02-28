# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PlanePartition(normal; [tol])

A method for partitioning objects into a family of hyperplanes defined by
a `normal` direction. Two points `pᵢ` and `pⱼ` belong to the same hyperplane
when `(pᵢ - pⱼ) ⋅ normal < tol`.
"""
struct PlanePartition{V<:Vec,ℒ<:Len} <: PointPredicatePartitionMethod
  normal::V
  tol::ℒ
  PlanePartition(normal::V, tol::ℒ) where {V<:Vec,ℒ<:Len} = new{V,float(ℒ)}(unormalize(normal), tol)
end

PlanePartition(normal::Vec, tol) = PlanePartition(normal, addunit(tol, u"m"))

PlanePartition(normal::Vec; tol=atol(eltype(normal))) = PlanePartition(normal, tol)

PlanePartition(normal::Tuple; kwargs...) = PlanePartition(Vec(normal); kwargs...)

(p::PlanePartition)(pᵢ, pⱼ) = abs(udot(pᵢ - pⱼ, p.normal)) < p.tol
