# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Tetrahedron(p1, p2, p3, p4)

A tetrahedron with points `p1`, `p2`, `p3`, `p4`.
"""
@polytope Tetrahedron 3 4

nvertices(::Type{<:Tetrahedron}) = 4

==(t1::Tetrahedron, t2::Tetrahedron) = t1.vertices == t2.vertices

Base.isapprox(t1::Tetrahedron, t2::Tetrahedron; kwargs...) =
  all(isapprox(v1, v2; kwargs...) for (v1, v2) in zip(t1.vertices, t2.vertices))

function measure(t::Tetrahedron)
  A, B, C, D = t.vertices
  abs((A - D) ⋅ ((B - D) × (C - D))) / 6
end

function boundary(t::Tetrahedron)
  indices = [(3, 2, 1), (4, 1, 2), (4, 3, 1), (4, 2, 3)]
  SimpleMesh(pointify(t), connect.(indices))
end

function (t::Tetrahedron)(u, v, w)
  z = (1 - u - v - w)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1) || (z < 0 || z > 1)
    throw(DomainError((u, v, w), "invalid barycentric coordinates for tetrahedron."))
  end
  v₁, v₂, v₃, v₄ = coordinates.(t.vertices)
  Point(v₁ * z + v₂ * u + v₃ * v + v₄ * w)
end
