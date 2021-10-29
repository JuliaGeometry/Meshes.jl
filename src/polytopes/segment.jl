# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment with end points `p1`, `p2`.
The segment can be called as `s(t)` with `t` between
`0` and `1` to interpolate linearly between its endpoints.

See also [`Line`](@ref).
"""
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{1,Dim,T}
  vertices::V
end

isconvex(::Type{<:Segment}) = true

nvertices(::Type{<:Segment}) = 2
nvertices(s::Segment) = nvertices(typeof(s))

Base.minimum(s::Segment) = s.vertices[1]
Base.maximum(s::Segment) = s.vertices[2]
Base.extrema(s::Segment) = s.vertices[1], s.vertices[2]
measure(s::Segment) = norm(s.vertices[2] - s.vertices[1])
Base.length(s::Segment) = measure(s)

function (s::Segment)(t)
  if t < 0 || t > 1
    throw(DomainError(t, "s(t) is not defined for t outside [0, 1]."))
  end
  a, b = s.vertices
  a + t * (b - a)
end

function collinear(p1::Point{2,T}, p2::Point{2,T}, p3::Point{2,T}) where {T}
  p1, p2, p3 = coordinates.((p1, p2, p3))
  # comparing the slopes p1 to p2, and p2 to p3
  (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p3[2] - p2[2]) * (p2[1] - p1[1]) == 0
end

function Base.in(p::Point{2,T}, s::Segment{2,T}) where {T}
  a, b = s.vertices
  # i)  collinearity
  arecollinear = collinear(a, p, b)
  # (ii) compare dot product
  if arecollinear
    kap = (b - a) ⋅ (p - a)
    kab = (b - a) ⋅ (b - a)
    return 0 ≤ kap ≤ kab
  else
    return false
  end
end

