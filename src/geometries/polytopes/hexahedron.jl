# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Hexahedron(p1, p2, ..., p8)

A hexahedron with points `p1`, `p2`, ..., `p8`.
"""
@polytope Hexahedron 3 8

nvertices(::Type{<:Hexahedron}) = 8

==(h₁::Hexahedron, h₂::Hexahedron) = h₁.vertices == h₂.vertices

Base.isapprox(h₁::Hexahedron, h₂::Hexahedron; atol=atol(lentype(h₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(h₁.vertices, h₂.vertices))

function (h::Hexahedron)(u, v, w)
  A1, A2, A4, A3, A5, A6, A8, A7 = to.(h.vertices)
  withcrs(
    h,
    (1 - u) * (1 - v) * (1 - w) * A1 +
    u * (1 - v) * (1 - w) * A2 +
    (1 - u) * v * (1 - w) * A3 +
    u * v * (1 - w) * A4 +
    (1 - u) * (1 - v) * w * A5 +
    u * (1 - v) * w * A6 +
    (1 - u) * v * w * A7 +
    u * v * w * A8
  )
end
