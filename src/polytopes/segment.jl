# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment with end points `p1`, `p2`.
The segment can be called as `s(t)` with `t` between
`0` and `1` to interpolate linearly between its endpoints.

See also [`Rope`](@ref), [`Ring`](@ref), [`Line`](@ref).
"""
struct Segment{Dim,T} <: Chain{Dim,T}
  vertices::NTuple{2,Point{Dim,T}}
  Segment{Dim,T}(vertices::NTuple{2,Point{Dim,T}}) where {Dim,T} = new(vertices)
end

function Segment(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T}
  if length(vertices) ≠ 2
    throw(ArgumentError("The number of vertices must be 2"))
  end
  Segment{Dim,T}((first(vertices), last(vertices)))
end

isconvex(::Type{<:Segment}) = true

isclosed(::Type{<:Segment}) = false

isparametrized(::Type{<:Segment}) = true

nvertices(::Type{<:Segment}) = 2

vertices(s::Segment) = collect(s.vertices)

Base.minimum(s::Segment) = s.vertices[1]

Base.maximum(s::Segment) = s.vertices[2]

Base.extrema(s::Segment) = s.vertices[1], s.vertices[2]

measure(s::Segment) = norm(s.vertices[2] - s.vertices[1])

boundary(s::Segment) = PointSet(vertices(s))

center(s::Segment) = s(0.5)

function Base.in(p::Point{Dim,T}, s::Segment{Dim,T}) where {Dim,T}
  # given collinear points (a, b, p), the point p intersects
  # segment ab if and only if vectors satisfy 0 ≤ ap ⋅ ab ≤ ||ab||²
  a, b = s.vertices
  ab, ap = b - a, p - a
  iscollinear(a, b, p) && zero(T) ≤ ab ⋅ ap ≤ ab ⋅ ab
end

function Base.isapprox(s1::Segment, s2::Segment)
  v1, v2 = s1.vertices, s2.vertices
  isapprox(v1[1], v2[1]) && isapprox(v1[2], v2[2])
end

function (s::Segment)(t)
  if t < 0 || t > 1
    throw(DomainError(t, "s(t) is not defined for t outside [0, 1]."))
  end
  a, b = s.vertices
  a + t * (b - a)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Segment{Dim,T}}) where {Dim,T} =
  Segment(rand(rng, Point{Dim,T}, 2))
