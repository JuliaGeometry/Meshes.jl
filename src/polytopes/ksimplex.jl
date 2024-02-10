"""
    Simplex(p1, p2, ..., pK, pK+1)

A K-simplex is a "simplical" polytope that occupies space in K dimensions using K+1 points,
and is embedded in some dimension `Dim`.
For example, a line segment in 3-D space would be Simplex{1, 3} as it occupies a single dimension,
is made up of two points, and lies in 3d space.

Note that we only allow constructors of the form `Simplex(p1,p2,...,pk+1)`, not `Simplex(Tuple(...))`.
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

"Generate normal vector to every facet of simplex in (K+1) dimensions."
function normal(splx::Simplex{K,Dim,T}) where {K, Dim, T<:Real}
    verts = vertices(splx)
    p0 = first(verts)
    extended_basis = [(p .- p0 for p in verts[2:end])... rand!(similar(p0.coords))]
    normal = qr(extended_basis).Q[:, end]
end

function measure(splx::Simplex{K,Dim,T}) where {K, Dim, T<:Real}
    # https://en.wikipedia.org/wiki/Cayley%E2%80%93Menger_determinant
    Ds_ = pairwise(SqEuclidean(), getfield.(vertices(splx), :coords))
    Ds = [Ds_ ones(size(Ds_, 1), 1);
          ones(1, size(Ds_, 2)) 0]
    factor = (-1)^(K+1)/(factorial(K)^2*2^K)
    return sqrt(factor * det(Ds))
end
