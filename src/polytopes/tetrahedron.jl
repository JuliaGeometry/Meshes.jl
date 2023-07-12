# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Tetrahedron(p1, p2, p3, p4)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`.
"""
Tetrahedron

@polytope Tetrahedron 3 4

issimplex(::Type{<:Tetrahedron}) = true

isconvex(::Type{<:Tetrahedron}) = true

isparametrized(::Type{<:Tetrahedron}) = true

nvertices(::Type{<:Tetrahedron}) = 4

vertices(t::Tetrahedron) = collect(t.vertices)

function measure(t::Tetrahedron)
  A, B, C, D = t.vertices
  abs((A - D) ⋅ ((B - D) × (C - D))) / 6
end

function boundary(t::Tetrahedron)
  indices = [(3, 2, 1), (4, 1, 2), (4, 3, 1), (4, 2, 3)]
  SimpleMesh(vertices(t), connect.(indices))
end

function (t::Tetrahedron)(u, v, w)
  z = (1 - u - v - w)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1) || (z < 0 || z > 1)
    throw(DomainError((u, v, w), "invalid barycentric coordinates for tetrahedron."))
  end
  v₁, v₂, v₃, v₄ = coordinates.(t.vertices)
  Point(v₁ * z + v₂ * u + v₃ * v + v₄ * w)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Tetrahedron{Dim,T}}) where {Dim,T} =
  Tetrahedron(rand(rng, Point{Dim,T}, 4))
