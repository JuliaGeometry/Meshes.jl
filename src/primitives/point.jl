# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

ncoords(::Type{<:CoordRefSystems.Geographic}) = 2
ncoords(::Type{<:CoordRefSystems.Projected}) = 2
ncoords(::Type{<:LatLonAlt}) = 3
ncoords(::Type{<:Cartesian{Datum,N}}) where {Datum,N} = N
ncoords(::Type{<:Polar}) = 2
ncoords(::Type{<:Cylindrical}) = 3
ncoords(::Type{<:Spherical}) = 3

asvec(coords::Cartesian) = Vec(coords.coords)
ascart(vec::Vec) = Cartesian(Tuple(vec))

struct Point{C<:CRS}
  coords::C
end

# convenience constructors
Point(coords...) = Point(Cartesian(coords...))

embeddim(::Type{Point{C}}) where {C<:CRS} = ncoords(C)
embeddim(::Type{Point{<:CoordRefSystems.Geographic}}) = 3

paramdim(::Type{<:Point}) = 0

center(p::Point) = p

==(A::Point, B::Point) = A.coords == B.coords

Base.isapprox(A::Point, B::Point; kwargs...) = isapprox(A.coords, B.coords; kwargs...)

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
-(A::Point{<:Cartesian}, B::Point{<:Cartesian}) = asvec(A.coords) - asvec(B.coords)

"""
    +(A::Point, v::Vec)
    +(v::Vec, A::Point)

Return the point at the end of the vector `v` placed
at a reference (or start) point `A`.
"""
+(A::Point{<:Cartesian}, v::Vec) = Point(ascart(asvec(A.coords) + v))
+(v::Vec, A::Point{<:Cartesian}) = A + v

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point{<:Cartesian}, v::Vec) = Point(ascart(asvec(A.coords) - v))
-(v::Vec, A::Point{<:Cartesian}) = A - v

"""
    ⪯(A::Point, B::Point)
    ⪰(A::Point, B::Point)
    ≺(A::Point, B::Point)
    ≻(A::Point, B::Point)

Generalized inequality for non-negative orthant Rⁿ₊.
"""
⪯(A::Point{<:Cartesian}, B::Point{<:Cartesian}) = all(≥(0), B - A)
⪰(A::Point{<:Cartesian}, B::Point{<:Cartesian}) = all(≥(0), A - B)
≺(A::Point{<:Cartesian}, B::Point{<:Cartesian}) = all(>(0), B - A)
≻(A::Point{<:Cartesian}, B::Point{<:Cartesian}) = all(>(0), A - B)

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
∠(A::P, B::P, C::P) where {P<:Point{<:Cartesian}} = ∠(A - B, C - B)

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Point{Cartesian{Datum,N,T}}}) where {Datum,N,T}
  tup = ntuple(Returns(rand(rng, T)), N)
  Point(Cartesian(tup))
end

# -----------
# IO METHODS
# -----------

# TODO: show
