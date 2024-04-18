# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

asvec(coords::Cartesian) = Vec(CoordRefSystems._coords(coords))
ascart(vec::Vec) = Cartesian(Tuple(vec))

struct Point{Dim,C<:CRS} <: Primitive{Dim}
  coords::C
  Point(coords::C) where {C<:CRS} = new{CoordRefSystems.ndims(coords),C}(coords)
end

# convenience constructors
Point(coords...) = Point(Cartesian(coords...))
Point(coords::Tuple) = Point(Cartesian(coords...))

paramdim(::Type{<:Point}) = 0

center(p::Point) = p

==(A::Point, B::Point) = A.coords == B.coords

Base.isapprox(A::Point, B::Point; atol=CoordRefSystems.tol(A.coords), kwargs...) =
  isapprox(A.coords, B.coords; atol, kwargs...)

"""
    coordinates(point)

Return the coordinates of the `point` with respect to the
canonical Euclidean basis.
"""
coordinates(A::Point) = A.coords

"""
    -(A::Point, B::Point)

Return the [`Vec`](@ref) associated with the direction
from point `B` to point `A`.
"""
-(A::Point{Dim,<:Cartesian}, B::Point{Dim,<:Cartesian}) where {Dim} = asvec(A.coords) - asvec(B.coords)

"""
    +(A::Point, v::Vec)
    +(v::Vec, A::Point)

Return the point at the end of the vector `v` placed
at a reference (or start) point `A`.
"""
+(A::Point{Dim,<:Cartesian}, v::Vec{Dim}) where {Dim} = Point(ascart(asvec(A.coords) + v))
+(v::Vec{Dim}, A::Point{Dim,<:Cartesian}) where {Dim} = A + v

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point{Dim,<:Cartesian}, v::Vec{Dim}) where {Dim} = Point(ascart(asvec(A.coords) - v))
-(v::Vec{Dim}, A::Point{Dim,<:Cartesian}) where {Dim} = A - v

"""
    ⪯(A::Point, B::Point)
    ⪰(A::Point, B::Point)
    ≺(A::Point, B::Point)
    ≻(A::Point, B::Point)

Generalized inequality for non-negative orthant Rⁿ₊.
"""
⪯(A::Point{Dim,<:Cartesian}, B::Point{Dim,<:Cartesian}) where {Dim} = all(≥(0u"m"), B - A)
⪰(A::Point{Dim,<:Cartesian}, B::Point{Dim,<:Cartesian}) where {Dim} = all(≥(0u"m"), A - B)
≺(A::Point{Dim,<:Cartesian}, B::Point{Dim,<:Cartesian}) where {Dim} = all(>(0u"m"), B - A)
≻(A::Point{Dim,<:Cartesian}, B::Point{Dim,<:Cartesian}) where {Dim} = all(>(0u"m"), A - B)

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

# TODO
# Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Point{Dim,T}}) where {Dim,T} = Point(rand(rng, Vec{Dim,T}))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, point::Point)
  if get(io, :compact, false)
    print(io, "(")
    _printfields(io, point.coords)
    print(io, ")")
  else
    print(io, "Point(")
    _printfields(io, point.coords)
    print(io, ")")
  end
end

_printfields(io, coords::CRS) = printfields(io, coords, compact=true)
_printfields(io, coords::CoordRefSystems.ShiftedCRS) = printfields(io, CoordRefSystems._coords(coords), compact=true)
_printfields(io, coords::Cartesian) =
  printfields(io, CoordRefSystems._coords(coords), CoordRefSystems._fnames(coords), compact=true)

function Base.show(io::IO, mime::MIME"text/plain", point::Point)
  print(io, "Point with ")
  show(io, mime, point.coords)
end
