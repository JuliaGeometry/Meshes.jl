# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Point(x₁, x₂, ..., xₙ)
    Point((x₁, x₂, ..., xₙ))

A point in `Dim`-dimensional space with coordinates in length units (default to meters).

The coordinates of the point are given with respect to the canonical
Euclidean basis, and integer coordinates are converted to float.

## Examples

```julia
# 2D points
Point(1.0, 2.0) # add default units
Point(1.0m, 2.0m) # double precision as expected
Point(1f0km, 2f0km) # single precision as expected
Point(1m, 2m) # integer is converted to float by design

# 3D points
Point(1.0, 2.0, 3.0) # add default units
Point(1.0m, 2.0m, 3.0m) # double precision as expected
Point(1f0km, 2f0km, 3f0km) # single precision as expected
Point(1m, 2m, 3m) # integer is converted to float by design
```

### Notes

Integer coordinates are not supported because most geometric processing
algorithms assume a continuous space. The conversion to float avoids
`InexactError` and other unexpected results.
"""
struct Point{M<:Manifold,C<:CRS} <: Primitive{M,C}
  coords::C
end

Point{M}(coords::C) where {M<:Manifold,C<:CRS} = Point{M,C}(coords)

Point(coords::CRS) = Point{_manifold(coords)}(coords)

# convenience constructor
Point(coords...) = Point(Cartesian(coords...))

# conversion
Base.convert(::Type{Point{M,CRSₜ}}, p::Point{M,CRSₛ}) where {M,CRSₜ,CRSₛ} = Point{M}(convert(CRSₜ, p.coords))
Base.convert(::Type{Point{M,CRS}}, p::Point{M,CRS}) where {M,CRS} = p

# promotion
function Base.promote(A::Point, B::Point)
  a, b = promote(A.coords, B.coords)
  Point(a), Point(b)
end

paramdim(::Type{<:Point}) = 0

function ==(A::Point, B::Point)
  A′, B′ = promote(A, B)
  to(A′) == to(B′)
end

function ==(A::Point{🌐}, B::Point{🌐})
  A′, B′ = promote(A, B)
  latlon₁ = convert(LatLon, A′.coords)
  latlon₂ = convert(LatLon, B′.coords)
  lat₁, lon₁ = latlon₁.lat, latlon₁.lon
  lat₂, lon₂ = latlon₂.lat, latlon₂.lon
  lat₁ == lat₂ && lon₁ == lon₂ || (abs(lon₁) == 180u"°" && lon₁ == -lon₂)
end

function Base.isless(A::Point, B::Point)
  A′, B′ = promote(A, B)
  isless(to(A′), to(B′))
end

function Base.isapprox(A::Point, B::Point; atol=atol(lentype(A)), kwargs...)
  A′, B′ = promote(A, B)
  isapprox(to(A′), to(B′); atol, kwargs...)
end

"""
    coords(point)

Return the coordinates of the `point`.
"""
coords(A::Point) = A.coords

"""
    to(point)

Return the vector from the origin to the `point`.
"""
to(A::Point) = Vec(CoordRefSystems.values(convert(Cartesian, A.coords)))

"""
    -(A::Point, B::Point)

Return the [`Vec`](@ref) associated with the direction
from point `B` to point `A`.
"""
-(A::Point, B::Point) = to(A) - to(B)

"""
    +(A::Point, v::Vec)
    +(v::Vec, A::Point)

Return the point at the end of the vector `v` placed
at a reference (or start) point `A`.
"""
+(A::Point, v::Vec) = withcrs(A, to(A) + v)
+(v::Vec, A::Point) = A + v

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point, v::Vec) = withcrs(A, to(A) - v)
-(v::Vec, A::Point) = A - v

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See <https://en.wikipedia.org/wiki/Angle>.

Uses the two-argument form of `atan` returning value in range [-π, π]
in 2D and [0, π] in 3D.
See <https://en.wikipedia.org/wiki/Atan2>.

## Examples

```julia
∠(Point(1,0), Point(0,0), Point(0,1)) == π/2
```
"""
∠(A::P, B::P, C::P) where {P<:Point} = ∠(A - B, C - B)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, point::Point)
  if get(io, :compact, false)
    print(io, "(")
  else
    print(io, "Point(")
  end
  values = CoordRefSystems.values(point.coords)
  names = CoordRefSystems.names(point.coords)
  printfields(io, values, names, singleline=true)
  print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", point::Point)
  print(io, "Point with ")
  show(io, mime, point.coords)
end

# -----------------
# HELPER FUNCTIONS
# -----------------

_manifold(coords::CRS) = 𝔼{CoordRefSystems.ndims(coords)}
_manifold(::LatLon) = 🌐
_manifold(::GeocentricLatLon) = 🌐
_manifold(::AuthalicLatLon) = 🌐
