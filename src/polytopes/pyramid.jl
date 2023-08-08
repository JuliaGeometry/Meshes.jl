# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pyramid(p1, p2, p3, p4, p5)

A pyramid with points `p1`, `p2`, `p3`, `p4`, `p5`.
"""
@polytope Pyramid 3 5

nvertices(::Type{<:Pyramid}) = 5

==(p1::Pyramid, p2::Pyramid) = p1.vertices == p2.vertices

Base.isapprox(p1::Pyramid, p2::Pyramid; kwargs...) =
  all(isapprox(v1, v2; kwargs...) for (v1, v2) in zip(p1.vertices, p2.vertices))

function boundary(p::Pyramid)
  indices = [(4, 3, 2, 1), (5, 1, 2), (5, 4, 1), (5, 3, 4), (5, 2, 3)]
  SimpleMesh(pointify(p), connect.(indices))
end
