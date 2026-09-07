# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented line segment from point `p1` to point `p2`.

See also [`Rope`](@ref), [`Ring`](@ref), [`Line`](@ref).
"""
struct Segment{M<:Manifold,C<:CRS,V<:Vec,ℒ<:Len} <: Chain{M,C}
  # input fields
  a::Point{M,C}
  b::Point{M,C}

  # state fields
  v::V
  l::ℒ
end

function Segment(a::Point, b::Point)
  a′, b′ = promote(a, b)
  v = _geodesicvec(a′, b′)
  l = _geodesiclen(a′, b′)
  Segment(a′, b′, v, l)
end

Segment(a::Tuple, b::Tuple) = Segment(Point(a), Point(b))

Segment(ab::AbstractVector) = Segment(ab...)

Segment(ab::Tuple) = Segment(ab...)

nvertices(::Type{<:Segment}) = 2

Base.minimum(s::Segment) = s.a

Base.maximum(s::Segment) = s.b

Base.extrema(s::Segment) = s.a, s.b

Base.length(s::Segment) = s.l

center(s::Segment{<:𝔼}) = s(1 // 2)

==(s₁::Segment, s₂::Segment) = s₁.a == s₂.a && s₁.b == s₂.b

Base.isapprox(s₁::Segment, s₂::Segment; atol=atol(lentype(s₁)), kwargs...) =
  isapprox(s₁.a, s₂.a; atol=atol, kwargs...) && isapprox(s₁.b, s₂.b; atol=atol, kwargs...)

(s::Segment{<:𝔼})(t) = s.a + (t * ustrip(s.l)) * s.v

(s::Segment{🌐})(t) = geodesicfwd(s.a, geodesicazimuth(s.a, s.v), t * s.l)

Base.reverse(s::Segment) = Segment(s.b, s.a)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, s::Segment)
  name = prettyname(s)
  ioctx = IOContext(io, :compact => true)
  print(io, "$name(")
  printfields(ioctx, s, (:a, :b), singleline=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::Segment)
  summary(io, s)
  printfields(io, s, (:a, :b))
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_geodesicvec(a::Point{𝔼{Dim}}, b::Point{𝔼{Dim}}) where {Dim} = unormalize(b - a)
_geodesicvec(a::Point{🌐}, b::Point{🌐}) = geodesictangent(a, geodesicbwd(a, b))
_geodesiclen(a::Point, b::Point) = GeodesicDistance()(a, b)
