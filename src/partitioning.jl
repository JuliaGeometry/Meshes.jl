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
    PredicatePartitionMethod

A method for partitioning domain/data objects with predicate functions.
"""
abstract type PredicatePartitionMethod <: PartitionMethod end

function partitioninds(rng::AbstractRNG, domain::Domain, method::PredicatePartitionMethod)
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
    SPredicatePartitionMethod

A method for partitioning domain/data objects with spatial predicate functions.
"""
abstract type SPredicatePartitionMethod <: PartitionMethod end

function partitioninds(rng::AbstractRNG, domain::Domain, method::SPredicatePartitionMethod)
  nelms = nelements(domain)
  subsets = Vector{Int}[]
  for i in randperm(rng, nelms)
    p = centroid(domain, i)
    x = to(p)
    inserted = false
    for subset in subsets
      q = centroid(domain, subset[1])
      y = to(q)
      if method(x, y)
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
include("partitioning/predicate.jl")
include("partitioning/spatialpredicate.jl")
include("partitioning/product.jl")
include("partitioning/hierarchical.jl")
