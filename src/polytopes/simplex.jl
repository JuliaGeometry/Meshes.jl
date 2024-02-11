"""
    Simplex(p1, p2, ..., pK, pK+1)

A K-simplex is a polytope with parametric dimension K and K+1 vertices.
It generalizes the [`Segment`](@ref), the  [`Triangle`](@ref) and the
[`Tetrahedron`](@ref) to higher dimensions (K > 3).

See also [`issimplex`](@ref).

Note that we only allow constructors of the form `Simplex(p1,p2,...,pk+1)`, not `Simplex(Tuple(...))`.
"""
struct Simplex{K,Dim,T,N} <: Polytope{K,Dim,T}
  vertices::NTuple{N,Point{Dim,T}}
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

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(::Type{<:Simplex{K}}) where {K} = K+1

function Base.isapprox(p₁::SimplexT, p₂::SimplexT; kwargs...) where {SimplexT<:Simplex}
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(vertices(p₁), vertices(p₂)))
end
