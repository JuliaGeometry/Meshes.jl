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
  if t < 0 || t > 1
    throw(DomainError(t, "w(u, v, w) is not defined for coordinates outside [0, 1]."))
  end
  a1, a2, a3, b1, b2, b3 = w.vertices
  a = Triangle(a1, a2, a3)
  b = Triangle(b1, b2, b3)
  Segment(a(u, v), b(u, v))(w)
end