# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
struct Pyramid{Dim,T} <: Polyhedron{Dim,T}
  vertices::NTuple{5,Point{Dim,T}}
  Pyramid{Dim,T}(vertices::NTuple{5,Point{Dim,T}}) where {Dim,T} = new(vertices)
end

function Pyramid{Dim,T}(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T}
  N = length(vertices)
  N == 5 || throw(ArgumentError("Invalid number of vertices for Pyramid. Expected 5, got $N."))
  v = ntuple(i -> @inbounds(vertices[i]), 5)
  Pyramid{Dim,T}(v)
end

Pyramid(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T} = Pyramid{Dim,T}(vertices)

nvertices(::Type{<:Pyramid}) = 5
nvertices(p::Pyramid) = nvertices(typeof(p))

vertices(p::Pyramid) = collect(p.vertices)

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(vertices(p), connect.(indices))
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Pyramid{Dim,T}}) where {Dim,T} =
  Pyramid(rand(rng, Point{Dim,T}, 5))
