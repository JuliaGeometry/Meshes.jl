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
  dirvec::V
  length::ℒ
end

function Geodesic(a::Point, b::Point)
  a′, b′ = promote(a, b)
  dirvec = _geodesicdirection(a′, b′)
  length = GeodesicDistance()(a′, b′)
  Geodesic(a′, b′, dirvec, length)
end

Geodesic(a::Tuple, b::Tuple) = Geodesic(Point(a), Point(b))

paramdim(::Type{<:Geodesic}) = 1

Base.minimum(g::Geodesic) = g.a

Base.maximum(g::Geodesic) = g.b

Base.extrema(g::Geodesic) = g.a, g.b

==(g₁::Geodesic, g₂::Geodesic) = g₁.a == g₂.a && g₁.b == g₂.b

Base.isapprox(g₁::Geodesic, g₂::Geodesic; atol=atol(lentype(g₁)), kwargs...) =
  isapprox(g₁.a, g₂.a; atol=atol, kwargs...) && isapprox(g₁.b, g₂.b; atol=atol, kwargs...)

(g::Geodesic{𝔼})(t) = g.a + (t * g.length) * g.dirvec

(g::Geodesic{🌐})(t) = geodesicfwd(g.a, geodesicazimuth(g.a, g.dirvec), t * g.length)

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

_geodesicdirection(a::Point{𝔼}, b::Point{𝔼}) = unormalize(b - a)
_geodesicdirection(a::Point{🌐}, b::Point{🌐}) = geodesictangent(a, geodesicbwd(a, b))
