# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PartitionMethod

A method for partitioning domain/data objects.
"""
abstract type PartitionMethod end

"""
    partition([rng], object, method)

Partition `object` with partition `method`.
Optionally, specify random number generator `rng`.
"""
partition(object, method::PartitionMethod) = partition(Random.default_rng(), object, method)

function partition(rng::AbstractRNG, object, method::PartitionMethod)
  subsets, metadata = partitioninds(rng, object, method)
  Partition(object, subsets, metadata)
end

"""
    partitioninds(rng, object, method)

Return subsets and metadata for the partition `method`
applied to the `object` with random number generator `rng`.
"""
function partitioninds end

"""
    IndexPredicatePartitionMethod

A method for partitioning a domain/data object with predicate functions
of the form `pred(i, j)` where `i` and `j` are the linear indices of the
elements of the object.
"""
abstract type IndexPredicatePartitionMethod <: PartitionMethod end

function partitioninds(rng::AbstractRNG, domain::Domain, method::IndexPredicatePartitionMethod)
  nelms = nelements(domain)
  subsets = Vector{Int}[]
  for i in randperm(rng, nelms)
    inserted = false
    for subset in subsets
      j = subset[1]
      if method(i, j)
        push!(subset, i)
        inserted = true
        break
      end
    end
    if !inserted
      push!(subsets, [i])
    end
  end

  subsets, Dict()
end

"""
    PointPredicatePartitionMethod

A method for partitioning a domain/data object with predicate functions
of the form `pred(pᵢ, pⱼ)` where `pᵢ` and `pⱼ` are the centroid points
of the the `i`-th and `j-th` elements of the object.
"""
abstract type PointPredicatePartitionMethod <: PartitionMethod end

function partitioninds(rng::AbstractRNG, domain::Domain, method::PointPredicatePartitionMethod)
  nelms = nelements(domain)
  subsets = Vector{Int}[]
  for i in randperm(rng, nelms)
    pᵢ = centroid(domain, i)
    inserted = false
    for subset in subsets
      pⱼ = centroid(domain, subset[1])
      if method(pᵢ, pⱼ)
        push!(subset, i)
        inserted = true
        break
      end
    end
    if !inserted
      push!(subsets, [i])
    end
  end

  subsets, Dict()
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("partitioning/uniform.jl")
include("partitioning/fraction.jl")
include("partitioning/block.jl")
include("partitioning/bisectpoint.jl")
include("partitioning/bisectfraction.jl")
include("partitioning/ball.jl")
include("partitioning/plane.jl")
include("partitioning/direction.jl")
include("partitioning/indexpred.jl")
include("partitioning/pointpred.jl")
include("partitioning/product.jl")
include("partitioning/hierarchical.jl")
