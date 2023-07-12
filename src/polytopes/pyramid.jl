# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
Pyramid

@polytope Pyramid 3 5

nvertices(::Type{<:Pyramid}) = 5

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(collect(p.vertices), connect.(indices))
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Pyramid{Dim,T}}) where {Dim,T} =
  Pyramid(rand(rng, Point{Dim,T}, 5))
