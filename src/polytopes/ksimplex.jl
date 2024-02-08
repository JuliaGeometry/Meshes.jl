"""
    Ksimplex(p1, p2, ..., pK, pK+1)

A K-simplex is a "simplical" polytope that occupies space in K dimensions using K+1 points,
and is embedded in some dimension `Dim`.
For example, a line segment in 3-D space would be Ksimplex{1, 3} as it occupies a single dimension,
is made up of two points, and lies in 3d space.
"""

@computed struct Ksimplex{K,Dim,T} <: Polytope{K,Dim,T}
    vertices::NTuple{K+1, Point{Dim, T}}
end

Ksimplex(vertices::Vararg{Tuple,N}) where {N} = Ksimplex(Point.(vertices))
Ksimplex(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ksimplex{N-1,Dim,T}(vertices)

Ksimplex{N}(vertices::Vararg{Tuple,N}) where {N} = Ksimplex(Point.(vertices))
Ksimplex{N}(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ksimplex{N-1,Dim,T}(vertices)
Ksimplex{N}(vertices::NTuple{N,Point{Dim,T}}) where {N,Dim,T} = Ksimplex{N-1,Dim,T}(vertices)

nvertices(::Type{<:Ksimplex{K}}) where {K} = K+1

function Base.isapprox(p₁::KSimplexT, p₂::KSimplexT; kwargs...) where {KSimplexT<:Ksimplex}
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))
end
