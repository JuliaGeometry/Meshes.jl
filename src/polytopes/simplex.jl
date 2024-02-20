# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Simplex(p1, p2, ..., pK, pK+1)

A K-simplex is a polytope with parametric dimension K and K+1 vertices.
It generalizes the [`Segment`](@ref), the  [`Triangle`](@ref) and the
[`Tetrahedron`](@ref) to higher dimensions (K > 3).

See also [`issimplex`](@ref).
"""
struct Simplex{K,Dim,T,N} <: Polytope{K,Dim,T}
  vertices::NTuple{N,Point{Dim,T}}
end

function Simplex(vertices::NTuple{N,Point{Dim,T}}) where {Dim,T,N}
  K = N - 1
  K ≤ Dim || throw(ArgumentError("rank (number of vertices minus one) must be less or equal to embedding dimension"))
  Simplex{K,Dim,T,N}(vertices)
end

Simplex(vertices::Point{Dim,T}...) where {Dim,T} = Simplex(vertices)
Simplex(vertices::Tuple...) = Simplex(Point.(vertices))

function Simplex{K}(vertices::NTuple{N,Point{Dim,T}}) where {K,Dim,T,N}
  N == K + 1 || throw(ArgumentError("number of vertices must be equal to rank K plus one"))
  K ≤ Dim || throw(ArgumentError("rank K must be less or equal to embedding dimension"))
  Simplex{K,Dim,T,N}(vertices)
end

Simplex{K}(vertices::Point{Dim,T}...) where {K,Dim,T} = Simplex{K}(vertices)
Simplex{K}(vertices::Tuple...) where {K} = Simplex{K}(Point.(vertices))

nvertices(::Type{<:Simplex{K}}) where {K} = K + 1

Base.isapprox(s₁::Simplex, s₂::Simplex; kwargs...) =
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(vertices(s₁), vertices(s₂)))
