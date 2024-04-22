# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PlanePartition(normal; tol=1e-6)

A method for partitioning spatial objects into a family of hyperplanes
defined by a `normal` direction. Two points `x` and `y` belong to the same
hyperplane when `(x - y) ⋅ normal < tol`.
"""
struct PlanePartition{V<:Vec} <: SPredicatePartitionMethod
  normal::V
  tol::Float64

  function PlanePartition{V}(normal, tol) where {V<:Vec}
    new(normalize(normal), tol)
  end
end

PlanePartition(normal::V; tol=1e-6) where {V<:Vec} = PlanePartition{V}(normal, tol)

PlanePartition(normal::Tuple; tol=1e-6) = PlanePartition(Vec(normal), tol=tol)

(p::PlanePartition)(x, y) = abs((x - y) ⋅ p.normal) < p.tol
