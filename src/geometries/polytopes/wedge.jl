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

function (w::Wedge)(u, v, w)
  ℒ = lentype(w)
  T = numtype(ℒ)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v, w), "w(u, v, w) is not defined for u, v, w outside [0, 1]³."))
  end
  a1, a2, a3, b1, b2, b3 = w.vertices
  a = Triangle(a1, a2, a3)
  b = Triangle(b1, b2, b3)
  s = Segment(a(T(u), T(v)), b(T(u), T(v)))
  s(T(w))
end