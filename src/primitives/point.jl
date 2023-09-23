# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Point(x₁, x₂, ..., xₙ)
    Point((x₁, x₂, ..., xₙ))
    Point{Dim,T}(x₁, x₂, ..., xₙ)
    Point{Dim,T}((x₁, x₂, ..., xₙ))

A point in `Dim`-dimensional space with coordinates of type `T`.

The coordinates of the point are given with respect to the canonical
Euclidean basis, and `Integer` coordinates are converted to `Float64`.

## Examples

```julia
# 2D points
A = Point(0.0, 1.0) # double precision as expected
B = Point(0f0, 1f0) # single precision as expected
C = Point(0, 0) # Integer is converted to Float64 by design
D = Point2(0, 1) # explicitly ask for double precision
E = Point2f(0, 1) # explicitly ask for single precision

# 3D points
F = Point(1.0, 2.0, 3.0) # double precision as expected
G = Point(1f0, 2f0, 3f0) # single precision as expected
H = Point(1, 2, 3) # Integer is converted to Float64 by design
I = Point3(1, 2, 3) # explicitly ask for double precision
J = Point3f(1, 2, 3) # explicitly ask for single precision
```

### Notes

- Type aliases are `Point1`, `Point2`, `Point3`, `Point1f`, `Point2f`, `Point3f`
- `Integer` coordinates are not supported because most geometric processing
  algorithms assume a continuous space. The conversion to `Float64` avoids
  `InexactError` and other unexpected results.
"""
struct Point{Dim,T} <: Primitive{Dim,T}
  coords::Vec{Dim,T}
  Point(coords::Vec{Dim,T}) where {Dim,T} = new{Dim,T}(coords)
end

# convenience constructors
Point{Dim,T}(coords...) where {Dim,T} = Point(Vec{Dim,T}(coords...))
Point(coords...) = Point(Vec(coords...))

# coordinate type conversions
Base.convert(::Type{Point{Dim,T}}, coords) where {Dim,T} = Point{Dim,T}(coords)
Base.convert(::Type{Point{Dim,T}}, p::Point) where {Dim,T} = Point{Dim,T}(p.coords)
Base.convert(::Type{Point}, coords) = Point{length(coords),eltype(coords)}(coords)

# type aliases for convenience
const Point1 = Point{1,Float64}
const Point2 = Point{2,Float64}
const Point3 = Point{3,Float64}
const Point1f = Point{1,Float32}
const Point2f = Point{2,Float32}
const Point3f = Point{3,Float32}

paramdim(::Type{Point{Dim,T}}) where {Dim,T} = 0

center(p::Point) = p

==(A::Point, B::Point) = A.coords == B.coords

Base.isapprox(A::Point{Dim,T}, B::Point{Dim,T}; atol=atol(T), kwargs...) where {Dim,T} =
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
-(A::Point, B::Point) = A.coords - B.coords

"""
    +(A::Point, v::Vec)
    +(v::Vec, A::Point)

Return the point at the end of the vector `v` placed
at a reference (or start) point `A`.
"""
+(A::Point, v::Vec) = Point(A.coords + v)
+(v::Vec, A::Point) = A + v

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point, v::Vec) = Point(A.coords - v)
-(v::Vec, A::Point) = A - v

"""
    ⪯(A::Point, B::Point)
    ⪰(A::Point, B::Point)
    ≺(A::Point, B::Point)
    ≻(A::Point, B::Point)

Generalized inequality for non-negative orthant Rⁿ₊.
"""
⪯(A::Point{Dim,T}, B::Point{Dim,T}) where {Dim,T} = all(≥(zero(T)), B - A)
⪰(A::Point{Dim,T}, B::Point{Dim,T}) where {Dim,T} = all(≥(zero(T)), A - B)
≺(A::Point{Dim,T}, B::Point{Dim,T}) where {Dim,T} = all(>(zero(T)), B - A)
≻(A::Point{Dim,T}, B::Point{Dim,T}) where {Dim,T} = all(>(zero(T)), A - B)

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

Uses the two-argument form of `atan` returning value in range [-π, π]
in 2D and [0, π] in 3D.
See https://en.wikipedia.org/wiki/Atan2.

## Examples

```julia
∠(Point(1,0), Point(0,0), Point(0,1)) == π/2
```
"""
∠(A::P, B::P, C::P) where {P<:Point{2}} = ∠(A - B, C - B)
∠(A::P, B::P, C::P) where {P<:Point{3}} = ∠(A - B, C - B)

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Point{Dim,T}}) where {Dim,T} = Point(rand(rng, Vec{Dim,T}))

function Base.show(io::IO, point::Point)
  if get(io, :compact, false)
    print(io, Tuple(point.coords))
  else
    print(io, "Point$(Tuple(point.coords))")
  end
end

Base.show(io::IO, ::MIME"text/plain", point::Point) = show(io, point)
