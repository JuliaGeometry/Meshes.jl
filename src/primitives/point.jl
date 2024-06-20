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
Point(1.0u"m", 2.0u"m") # double precision as expected
Point(1f0u"km", 2f0u"km") # single precision as expected
Point(1u"m", 2u"m") # integer is converted to float by design

# 3D points
Point(1.0, 2.0, 3.0) # add default units
Point(1.0u"m", 2.0u"m", 3.0u"m") # double precision as expected
Point(1f0u"km", 2f0u"km", 3f0u"km") # single precision as expected
Point(1u"m", 2u"m", 3u"m") # integer is converted to float by design
```

### Notes

- Integer coordinates are not supported because most geometric processing
  algorithms assume a continuous space. The conversion to float avoids
  `InexactError` and other unexpected results.
"""
struct Point{Dim,C<:CRS} <: Primitive{Dim,C}
  coords::C
  Point(coords::C) where {C<:CRS} = new{CoordRefSystems.ndims(coords),C}(coords)
end

# convenience constructor
Point(coords...) = Point(Cartesian(coords...))

# conversions
Base.convert(::Type{Point{Dim,CRSₜ}}, p::Point{Dim,CRSₛ}) where {Dim,CRSₜ,CRSₛ} = Point(convert(CRSₜ, p.coords))
Base.convert(::Type{Point{Dim,CRS}}, p::Point{Dim,CRS}) where {Dim,CRS} = p

paramdim(::Type{<:Point}) = 0

center(p::Point) = p

==(A::Point, B::Point) = to(A) == to(B)

Base.isapprox(A::Point, B::Point; atol=atol(lentype(A)), kwargs...) = isapprox(to(A), to(B); atol, kwargs...)

"""
    coords(point)

Return the coordinates of the `point`.
"""
coords(A::Point) = A.coords

"""
    to(point)

Return the vector from the origin to the `point`.
"""
to(A::Point) = Vec(CoordRefSystems.cvalues(convert(Cartesian, A.coords)))

"""
    -(A::Point, B::Point)

Return the [`Vec`](@ref) associated with the direction
from point `B` to point `A`.
"""
-(A::Point{Dim}, B::Point{Dim}) where {Dim} = to(A) - to(B)

"""
    +(A::Point, v::Vec)
    +(v::Vec, A::Point)

Return the point at the end of the vector `v` placed
at a reference (or start) point `A`.
"""
+(A::Point{Dim}, v::Vec{Dim}) where {Dim} = withdatum(A, to(A) + v)
+(v::Vec{Dim}, A::Point{Dim}) where {Dim} = A + v

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point{Dim}, v::Vec{Dim}) where {Dim} = withdatum(A, to(A) - v)
-(v::Vec{Dim}, A::Point{Dim}) where {Dim} = A - v

"""
    ⪯(A::Point, B::Point)
    ⪰(A::Point, B::Point)
    ≺(A::Point, B::Point)
    ≻(A::Point, B::Point)

Generalized inequality for non-negative orthant Rⁿ₊.
"""
⪯(A::Point{Dim}, B::Point{Dim}) where {Dim} = all(x -> x ≥ zero(x), B - A)
⪰(A::Point{Dim}, B::Point{Dim}) where {Dim} = all(x -> x ≥ zero(x), A - B)
≺(A::Point{Dim}, B::Point{Dim}) where {Dim} = all(x -> x > zero(x), B - A)
≻(A::Point{Dim}, B::Point{Dim}) where {Dim} = all(x -> x > zero(x), A - B)

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
∠(A::P, B::P, C::P) where {P<:Point{2}} = ∠(A - B, C - B)
∠(A::P, B::P, C::P) where {P<:Point{3}} = ∠(A - B, C - B)

Random.rand(rng::Random.AbstractRNG, ::Type{Point{Dim}}) where {Dim} = Point(rand(rng, Cartesian{NoDatum,Dim}))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, point::Point)
  if get(io, :compact, false)
    print(io, "(")
  else
    print(io, "Point(")
  end
  cvalues = CoordRefSystems.cvalues(point.coords)
  cnames = CoordRefSystems.cnames(point.coords)
  printfields(io, cvalues, cnames, singleline=true)
  print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", point::Point)
  print(io, "Point with ")
  show(io, mime, point.coords)
end
