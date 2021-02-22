# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PartitionMethod

A method for partitioning domain/data objects.
"""
abstract type PartitionMethod end

"""
    partition(object, method)

Partition `object` with partition `method`.
"""
function partition end

"""
    PredicatePartitionMethod

A method for partitioning domain/data objects with predicate functions.
"""
abstract type PredicatePartitionMethod <: PartitionMethod end

function partition(object, method::PredicatePartitionMethod)
  subsets = Vector{Int}[]
  for i in randperm(nelements(object))
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

  Partition(object, subsets)
end

"""
    SPredicatePartitionMethod

A method for partitioning domain/data objects with spatial predicate functions.
"""
abstract type SPredicatePartitionMethod <: PartitionMethod end

function partition(object, method::SPredicatePartitionMethod)
  Dim = embeddim(object)
  T = coordtype(object)

  # pre-allocate memory for coordinates
  x = MVector{Dim,T}(undef)
  y = MVector{Dim,T}(undef)

  subsets = Vector{Int}[]
  for i in randperm(nelements(object))
    coordinates!(x, object, i)
    inserted = false
    for subset in subsets
      coordinates!(y, object, subset[1])
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

  Partition(object, subsets)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("partitioning/random.jl")
include("partitioning/fraction.jl")
include("partitioning/block.jl")
include("partitioning/bisect_point.jl")
include("partitioning/bisect_fraction.jl")
include("partitioning/ball.jl")
include("partitioning/plane.jl")
include("partitioning/direction.jl")
include("partitioning/predicate.jl")
# include("partitioning/spatial_predicate.jl")
# include("partitioning/product.jl")
# include("partitioning/hierarchical.jl")
