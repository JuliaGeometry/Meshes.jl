# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type AbstractFace{N,T} end

struct NgonFace{N,T} <: AbstractFace{N,T}
    data::NTuple{N,T}
end

@propagate_inbounds Base.getindex(x::NgonFace, i::Integer) = x.data[i]
@propagate_inbounds Base.iterate(x::NgonFace) = iterate(x.data)
@propagate_inbounds Base.iterate(x::NgonFace, i) = iterate(x.data, i)
Base.length(::NgonFace{N,T}) where {N,T} = N

Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F

const LineFace = NgonFace{2,Int}
const TriangleFace = NgonFace{3,Int}
const QuadFace = NgonFace{4,Int}
const TetrahedronFace = NgonFace{4,Int}
