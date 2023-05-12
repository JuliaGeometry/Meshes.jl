# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BlockSampling(sides)

A method for sampling objects that are `sides` apart using a
[`BlockPartition`](@ref).

    BlockSampling(side₁, side₂, ..., sideₙ)

Alternatively, specify the sides `side₁`, `side₂`, ..., `sideₙ`.
"""
struct BlockSampling{S} <: DiscreteSamplingMethod
  sides::S
end

BlockSampling(sides...) = BlockSampling(sides)

function sample(::AbstractRNG, object, method::BlockSampling)
  Π = partition(object, BlockPartition(method.sides))
  view(object, first.(indices(Π)))
end
