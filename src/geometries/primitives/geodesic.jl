# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geodesic(a, b)

A geodesic from point `a` to `b`.

See also [`Segment`](@ref), [`Line`](@ref).
"""
struct Geodesic{M<:Manifold,C<:CRS,V<:Vec,ℒ<:Len} <: Primitive{M,C}
  # input fields
  a::Point{M,C}
  b::Point{M,C}

  # state fields
  v::V
  l::ℒ
end

function Geodesic(a::Point, b::Point)
  a′, b′ = promote(a, b)
  v = _geodesicdir(a′, b′)
  l = _geodesiclen(a′, b′)
  Geodesic(a′, b′, v, l)
end

Geodesic(a::Tuple, b::Tuple) = Geodesic(Point(a), Point(b))

paramdim(::Type{<:Geodesic}) = 1

Base.minimum(g::Geodesic) = g.a

Base.maximum(g::Geodesic) = g.b

Base.extrema(g::Geodesic) = g.a, g.b

Base.length(g::Geodesic) = g.l

==(g₁::Geodesic, g₂::Geodesic) = g₁.a == g₂.a && g₁.b == g₂.b

Base.isapprox(g₁::Geodesic, g₂::Geodesic; atol=atol(lentype(g₁)), kwargs...) =
  isapprox(g₁.a, g₂.a; atol=atol, kwargs...) && isapprox(g₁.b, g₂.b; atol=atol, kwargs...)

(g::Geodesic{<:𝔼})(t) = g.a + (t * g.l) * g.v

(g::Geodesic{🌐})(t) = geodesicfwd(g.a, geodesicazimuth(g.a, g.v), t * g.l)

Base.reverse(g::Geodesic) = Geodesic(g.b, g.a)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, g::Geodesic)
  name = prettyname(g)
  ioctx = IOContext(io, :compact => true)
  print(io, "$name(")
  printfields(ioctx, g, (:a, :b), singleline=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", g::Geodesic)
  summary(io, g)
  printfields(io, g, (:a, :b))
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_geodesicdir(a::Point{𝔼{Dim}}, b::Point{𝔼{Dim}}) where {Dim} = unormalize(b - a)
_geodesicdir(a::Point{🌐}, b::Point{🌐}) = geodesictangent(a, geodesicbwd(a, b))
_geodesiclen(a::Point, b::Point) = GeodesicDistance()(a, b)
