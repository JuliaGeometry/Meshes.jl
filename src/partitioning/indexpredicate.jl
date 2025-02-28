# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IndexPredicatePartition(pred)

A method for partitioning objects with a given `pred`icate function.
Two linear indices `i` and `j` are part of the same subset whenever
`pred(i, j) == true`
"""
struct IndexPredicatePartition <: IndexPredicatePartitionMethod
  pred::Function
end

(p::IndexPredicatePartition)(i, j) = p.pred(i, j)
