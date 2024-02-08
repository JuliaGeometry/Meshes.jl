"""
    Ksimplex(p1, p2, ..., pK, pK+1)

A K-simplex is a "simplical" polytope that occupies space in K dimensions using K+1 points,
and is embedded in some dimension `Dim`.
For example, a line segment in 3-D space would be Ksimplex{1, 3} as it occupies a single dimension,
is made up of two points, and lies in 3d space.
"""

struct Ksimplex{K,Dim,T,K_} <: Polytope{K,Dim,T}
    vertices::NTuple{K_, Point{Dim, T}}
end

Ksimplex(vertices::Vararg{NTuple{Dim, T},K_}) where {K_, Dim, T} = Ksimplex{K_-1, Dim, T, K_}(Point.(vertices))
# Ksimplex(NTuple{K_, Point{Dim, T}}) where {K_, Dim, T} = Ksimplex{K_-1, Dim, T, K_}(vertices)
Ksimplex(vertices::Vararg{Point{Dim,T},K_}) where {K_,Dim,T} = Ksimplex{K_-1,Dim,T,K_}(vertices)

Ksimplex{K}(vertices::Vararg{Tuple,K_}) where {K,K_} = Ksimplex(Point.(vertices))
Ksimplex{K}(vertices::Vararg{Point{Dim,T},K_}) where {K,Dim,T,K_} = Ksimplex{K_-1,Dim,T,K_}(vertices)
Ksimplex{K}(vertices::NTuple{K_,Point{Dim,T}}) where {K,Dim,T,K_} = Ksimplex{K_-1,Dim,T,K_}(vertices)

Ksimplex{K,Dim}(vertices::Vararg{Tuple,K_}) where {K,Dim,K_} = Ksimplex(Point.(vertices))
Ksimplex{K,Dim}(vertices::Vararg{Point{Dim,T},K_}) where {K,Dim,T,K_} = Ksimplex{K_-1,Dim,T,K_}(vertices)
Ksimplex{K,Dim}(vertices::NTuple{K_,Point{Dim,T}}) where {K,Dim,T,K_} = Ksimplex{K_-1,Dim,T,K_}(vertices)

nvertices(::Type{<:Ksimplex{K}}) where {K} = K+1

function Base.isapprox(p₁::KSimplexT, p₂::KSimplexT; kwargs...) where {KSimplexT<:Ksimplex}
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))
end
