# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PlanePartition(normal; tol=1e-6)

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

PlanePartition(normal::Vec; tol=1e-6u"m") = PlanePartition(normal, tol)

PlanePartition(normal::Tuple; tol=1e-6u"m") = PlanePartition(Vec(normal), tol)

(p::PlanePartition)(x, y) = abs((x - y) ⋅ p.normal) < p.tol^2
