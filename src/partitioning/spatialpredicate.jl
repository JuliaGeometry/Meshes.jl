# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SpatialPredicatePartition(predicate)

A method for partitioning spatial objects with a given spatial
`predicate` function. Two coordinates `x` and `y` are part of
the same subset whenever `predicate(x, y) == true`.
"""
struct SpatialPredicatePartition <: SPredicatePartitionMethod
  pred::Function
end

(p::SpatialPredicatePartition)(x, y) = p.pred(x, y)
