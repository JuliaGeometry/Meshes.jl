# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Wedge(p1, p2, p3, p4, p5, p6)

A Wedge with points `p1`, `p2`, `p3`, `p4`, `p5`, `p6`.
"""
@polytope Wedge 3 6

nvertices(::Type{<:Wedge}) = 6

==(t₁::Wedge, t₂::Wedge) = t₁.vertices == t₂.vertices

Base.isapprox(t₁::Wedge, t₂::Wedge; atol=atol(lentype(t₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(t₁.vertices, t₂.vertices))

function (wedge::Wedge)(u, v, w)
  ℒ = lentype(wedge)
  T = numtype(ℒ)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v, w), "wedge(u, v, w) is not defined for u, v, w outside [0, 1]³."))
  end
  a1, a2, a3, b1, b2, b3 = vertices(wedge)
  a = Quadrangle(a1, b1, b2, a2)
  b = Quadrangle(a1, b1, b3, a3)
  uv = T(u), T(v)
  s = Segment(a(uv...), b(uv...))
  s(T(w))
end
