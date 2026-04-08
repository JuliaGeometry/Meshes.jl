# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Wedge(p1, p2, p3, p4, p5, p6)

A Wedge with points `p1`, `p2`, `p3`, `p4`, `p5`, `p6`.
"""
@polytope Wedge 3 6

nvertices(::Type{<:Wedge}) = 6

==(w₁::Wedge, w₂::Wedge) = w₁.vertices == w₂.vertices

Base.isapprox(w₁::Wedge, w₂::Wedge; atol=atol(lentype(w₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(w₁.vertices, w₂.vertices))

function (wedge::Wedge)(u, v, w)
  a₁, a₂, a₃, b₁, b₂, b₃ = wedge.vertices
  a = Quadrangle(a₁, b₁, b₂, a₂)(u, v)
  b = Quadrangle(a₁, b₁, b₃, a₃)(u, v)
  Segment(promote(a, b)...)(w)
end
