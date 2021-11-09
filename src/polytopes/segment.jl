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

function Base.in(p::Point{Dim,T}, s::Segment{Dim,T}) where {Dim,T}
  a, b = s.vertices
  if Dim == 2
    _0 = zero(T)
  elseif Dim == 3
    _0 = Vec{3,T}(0, 0, 0)
  else
    throw(ErrorException("not implemented"))
  end
  # given collinear points (a, b, p), the point p intersects
  # segment ab if and only if vectors satisfy 0 ≤ ap ⋅ ab ≤ ||ab||²
  iscollinear = isapprox((b - a) × (p - a), _0, atol=atol(T)^2)
  iscollinear && zero(T) ≤ (b - a) ⋅ (p - a) ≤ (b - a) ⋅ (b - a)
end
