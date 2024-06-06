# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Tetrahedron(p1, p2, p3, p4)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`.
"""
@polytope Tetrahedron 3 4

nvertices(::Type{<:Tetrahedron}) = 4

Base.isapprox(t₁::Tetrahedron, t₂::Tetrahedron; kwargs...) =
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(t₁.vertices, t₂.vertices))

function (t::Tetrahedron)(u, v, w)
  z = (1 - u - v - w)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1) || (z < 0 || z > 1)
    throw(DomainError((u, v, w), "invalid barycentric coordinates for tetrahedron."))
  end
  v₁, v₂, v₃, v₄ = to.(t.vertices)
  withdatum(t, v₁ * z + v₂ * u + v₃ * v + v₄ * w)
end

Random.rand(rng::Random.AbstractRNG, ::Type{Tetrahedron{Dim}}) where {Dim} =
  Tetrahedron(ntuple(i -> rand(rng, Point{Dim}), 4))
