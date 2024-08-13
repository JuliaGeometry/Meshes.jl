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
