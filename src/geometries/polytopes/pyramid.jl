# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.

The first four points form the base and the last is the is the apex.
"""
@polytope Pyramid 3 5

nvertices(::Type{<:Pyramid}) = 5

function base(p::Pyramid)
  a, b, c, d, o = p.vertices
  Quadrangle(a, b, c, d)
end

==(p₁::Pyramid, p₂::Pyramid) = p₁.vertices == p₂.vertices

Base.isapprox(p₁::Pyramid, p₂::Pyramid; atol=atol(lentype(p₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))

function (pyramid::Pyramid)(u, v, w)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v, w), "pyramid(u, v, w) is not defined for u, v, w outside [0, 1]³."))
  end
  o = last(pyramid.vertices)
  q = base(pyramid)
  s = Segment(q(u, v), o)
  s(w)
end
