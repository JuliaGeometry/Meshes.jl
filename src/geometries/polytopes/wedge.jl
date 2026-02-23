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
  a₁, a₂, a₃, b₁, b₂, b₃ = vertices(wedge)
  q₁ = Quadrangle(a₁, b₁, b₂, a₂)
  q₂ = Quadrangle(a₁, b₁, b₃, a₃)
  s = Segment(q₁(u, v), q₂(u, v))
  s(w)
end
