"""
    Simplex(p1, p2, ..., pK, pK+1)

A K-simplex is a "simplical" polytope that occupies space in K dimensions using K+1 points,
and is embedded in some dimension `Dim`.
For example, a line segment in 3-D space would be Simplex{1, 3} as it occupies a single dimension,
is made up of two points, and lies in 3d space.
"""

struct Simplex{K,Dim,T,K_} <: Polytope{K,Dim,T}
    vertices::NTuple{K_, Point{Dim, T}}
end

Simplex(vertices::Vararg{NTuple{Dim, T},K_}) where {K_, Dim, T} = Simplex{K_-1, Dim, T, K_}(Point.(vertices))
# Simplex(NTuple{K_, Point{Dim, T}}) where {K_, Dim, T} = Simplex{K_-1, Dim, T, K_}(vertices)
Simplex(vertices::Vararg{Point{Dim,T},K_}) where {K_,Dim,T} = Simplex{K_-1,Dim,T,K_}(vertices)

Simplex{K}(vertices::Vararg{Tuple,K_}) where {K,K_} = Simplex(Point.(vertices))
Simplex{K}(vertices::Vararg{Point{Dim,T},K_}) where {K,Dim,T,K_} = Simplex{K_-1,Dim,T,K_}(vertices)
Simplex{K}(vertices::NTuple{K_,Point{Dim,T}}) where {K,Dim,T,K_} = Simplex{K_-1,Dim,T,K_}(vertices)

Simplex{K,Dim}(vertices::Vararg{Tuple,K_}) where {K,Dim,K_} = Simplex(Point.(vertices))
Simplex{K,Dim}(vertices::Vararg{Point{Dim,T},K_}) where {K,Dim,T,K_} = Simplex{K_-1,Dim,T,K_}(vertices)
Simplex{K,Dim}(vertices::NTuple{K_,Point{Dim,T}}) where {K,Dim,T,K_} = Simplex{K_-1,Dim,T,K_}(vertices)

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
