# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
@polytope Pyramid 3 5

nvertices(::Type{<:Pyramid}) = 5

Base.isapprox(p₁::Pyramid, p₂::Pyramid; kwargs...) =
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))

Random.rand(rng::Random.AbstractRNG, ::Type{Pyramid{Dim}}) where {Dim} = Pyramid(ntuple(i -> rand(rng, Point{Dim}), 5))
