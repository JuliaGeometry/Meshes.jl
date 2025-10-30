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
struct BlockSampling{Dim,ℒ<:Len} <: DiscreteSamplingMethod
  sides::NTuple{Dim,ℒ}
  BlockSampling(sides::NTuple{Dim,ℒ}) where {Dim,ℒ<:Len} = new{Dim,float(ℒ)}(sides)
end

BlockSampling(sides::NTuple{Dim,Len}) where {Dim} = BlockSampling(promote(sides...))

BlockSampling(sides::Tuple) = BlockSampling(aslen.(sides))

BlockSampling(sides...) = BlockSampling(sides)

function sampleinds(::AbstractRNG, d::Domain, method::BlockSampling)
  Π = partition(d, BlockPartition(method.sides))
  first.(indices(Π))
end
