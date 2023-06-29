# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Hexahedron(p1, p2, ..., p8)

A hexahedron with points `p1`, `p2`, ..., `p8`.
"""
struct Hexahedron{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  vertices::V
end

isperiodic(::Type{<:Hexahedron}) = (false, false, false)

isparametrized(::Type{<:Hexahedron}) = true

nvertices(::Type{<:Hexahedron}) = 8
nvertices(h::Hexahedron) = nvertices(typeof(h))

function boundary(h::Hexahedron)
  indices = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(h.vertices, connect.(indices))
end

function (h::Hexahedron)(u, v, w)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v, w), "h(u, v, w) is not defined for u, v, w outside [0, 1]³."))
  end
  A1, A2, A4, A3, A5, A6, A8, A7 = coordinates.(h.vertices)
  Point(
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

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Hexahedron{Dim,T}}) where {Dim,T} =
  Hexahedron(rand(rng, Point{Dim,T}, 8))
