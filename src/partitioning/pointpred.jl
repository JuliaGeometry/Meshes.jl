# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointPredicatePartition(pred)

A method for partitioning objects with a given point `pred`icate function.
Two points `pᵢ` and `pⱼ` are part of the same subset whenever `pred(pᵢ, pⱼ) == true`.
"""
struct PointPredicatePartition <: PointPredicatePartitionMethod
  pred::Function
end

(p::PointPredicatePartition)(pᵢ, pⱼ) = p.pred(pᵢ, pⱼ)
