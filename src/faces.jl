# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type AbstractFace{N,T} end
abstract type AbstractSimplexFace{N,T} <: AbstractFace{N,T} end
abstract type AbstractNgonFace{N,T} <: AbstractFace{N,T} end

struct SimplexFace{N,T} <: AbstractSimplexFace{N,T}
    data::NTuple{N,T}
end

@propagate_inbounds Base.getindex(x::SimplexFace, i::Integer) = x.data[i]
@propagate_inbounds Base.iterate(x::SimplexFace) = iterate(x.data)
@propagate_inbounds Base.iterate(x::SimplexFace, i) = iterate(x.data, i)
Base.length(::SimplexFace{N,T}) where {N,T} = N

const TetrahedronFace{T} = SimplexFace{4,T}
Face(::Type{<:SimplexFace{N}}, ::Type{T}) where {N,T} = SimplexFace{N,T}

struct NgonFace{N,T} <: AbstractNgonFace{N,T}
    data::NTuple{N,T}
end

@propagate_inbounds Base.getindex(x::NgonFace, i::Integer) = x.data[i]
@propagate_inbounds Base.iterate(x::NgonFace) = iterate(x.data)
@propagate_inbounds Base.iterate(x::NgonFace, i) = iterate(x.data, i)
Base.length(::NgonFace{N,T}) where {N,T} = N

const LineFace = NgonFace{2,Int}
const TriangleFace = NgonFace{3,Int}
const QuadFace = NgonFace{4,Int}
Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F
