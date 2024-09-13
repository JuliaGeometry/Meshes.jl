# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PlanePartition(normal; [tol])

A method for partitioning spatial objects into a family of hyperplanes
defined by a `normal` direction. Two points `x` and `y` belong to the same
hyperplane when `(x - y) ⋅ normal < tol`.
"""
struct PlanePartition{V<:Vec,ℒ<:Len} <: SPredicatePartitionMethod
  normal::V
  tol::ℒ
  PlanePartition(normal::V, tol::ℒ) where {V<:Vec,ℒ<:Len} = new{V,float(ℒ)}(unormalize(normal), tol)
end

PlanePartition(normal::Vec, tol) = PlanePartition(normal, addunit(tol, u"m"))

PlanePartition(normal::Vec; tol=atol(eltype(normal))) = PlanePartition(normal, tol)

PlanePartition(normal::Tuple; kwargs...) = PlanePartition(Vec(normal); kwargs...)

function (p::PlanePartition)(x, y)
  v = x - y
  n = withunit.(p.normal, unit(eltype(v)))
  abs(udot(v, n)) < p.tol
end
