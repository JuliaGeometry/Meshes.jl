# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BlockSampling(sides; neighbors=false)
    BlockSampling(side₁, side₂, ...; neighbors=false)

A method for sampling objects that are `sides` apart
using a [`BlockPartition`](@ref).
"""
struct BlockSampling{Dim,T} <: DiscreteSamplingMethod
  sides::SVector{Dim,T}
end

BlockSampling(sides::NTuple) = BlockSampling(SVector(sides))
BlockSampling(sides::Vararg) = BlockSampling(SVector(sides))

function sample(::AbstractRNG, object, method::BlockSampling)
  Π = partition(object, BlockPartition(method.sides))
  view(object, first.(indices(Π)))
end