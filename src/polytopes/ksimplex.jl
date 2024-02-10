"""
    Simplex(p1, p2, ..., pK, pK+1)

A K-simplex is a "simplical" polytope that occupies space in K dimensions using K+1 points,
and is embedded in some dimension `Dim`.
For example, a line segment in 3-D space would be Simplex{1, 3} as it occupies a single dimension,
is made up of two points, and lies in 3d space.
"""
struct Simplex{K,Dim,T,K_} <: Polytope{K,Dim,T}
    # Notice that this class has no default constructor (on purpose), since `K` can't be determined.
    # Instead we provide different constructors that also include some dimensionality checking.
    vertices::NTuple{K_, Point{Dim, T}}
end

Simplex(vertices::Vararg{Point{Dim,T},K_})           where {  Dim,T,K_   } = let
    (K_-1)<=Dim  || throw(ArgumentError("(Rank K)==(num vertices - 1) must be less or equal to embedding dimension Dim."))
    Simplex{K_-1,Dim,T,K_}(vertices)
end
Simplex{K}(vertices::Vararg{Point{Dim,T},K_})        where {K,Dim,T,K_   } = let
    K+1 == K_ || throw(ArgumentError("Num vertices must be rank K plus one."))
    K <= Dim || throw(ArgumentError("Rank K must be less or equal to embedding dimension Dim."))
    Simplex{K_-1,Dim,T,K_}(vertices)
end
Simplex{K,Dim}(vertices::Vararg{Point{Dim,T},K_})    where {K,Dim,T,K_   } = let
    K+1 == K_ || throw(ArgumentError("Num vertices must be rank K plus one."))
    K<=Dim || throw(ArgumentError("Rank K must be less or equal to embedding dimension Dim."))
    Simplex{K_-1,Dim,T,K_}(vertices)
end
Simplex{K,Dim,T}(vertices::Vararg{Point{Dim,T′},K_}) where {K,Dim,T,K_,T′} = let
    K+1 == K_ || throw(ArgumentError("Num vertices must be rank K plus one"))
    K<=Dim || throw(ArgumentError("Rank K must be less or equal to embedding dimension Dim."))
    Simplex{K_-1,Dim,T,K_}(vertices)
end

nvertices(::Type{<:Simplex{K}}) where {K} = K+1

function Base.isapprox(p₁::SimplexT, p₂::SimplexT; kwargs...) where {SimplexT<:Simplex}
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))
end
