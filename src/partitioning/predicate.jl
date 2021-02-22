# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PredicatePartition(pred)

A method for partitioning spatial objects with a given `predicate`.
Two locations `i` and `j` are part of the same subset whenever
`predicate(i, j) == true`
"""
struct PredicatePartition <: PredicatePartitionMethod
  pred::Function
end

(p::PredicatePartition)(i, j) = p.pred(i, j)
