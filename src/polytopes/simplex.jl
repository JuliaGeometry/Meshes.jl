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
  K ≤ Dim || throw(ArgumentError("Simplex rank (number of vertices - 1) must be less or equal to embedding dimension."))
  Simplex{K,Dim,T,N}(vertices)
end

Simplex(vertices::Point{Dim,T}...) where {Dim,T} = Simplex(vertices)
Simplex(vertices::Tuple...) = Simplex(Point.(vertices))

# Parametric constructor required for [`materialize`](@ref).
Simplex{K}(vertices::NTuple{N}) where {K,N} = let
  N == K+1 || throw(ArgumentError("Number of vertices must be rank K plus one."))
  Simplex(vertices...)
end
Simplex{K}(vertices::Vararg{Point,N}) where {K,N} = let
  N == K+1 || throw(ArgumentError("Number of vertices must be rank K plus one."))
  Simplex(vertices...)
end

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(::Type{<:Simplex{K}}) where {K} = K+1

function Base.isapprox(p₁::SimplexT, p₂::SimplexT; kwargs...) where {SimplexT<:Simplex}
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(vertices(p₁), vertices(p₂)))
end
