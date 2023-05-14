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
struct Point{Dim,T}
  coords::Vec{Dim,T}
  Point{Dim,T}(coords::Vec{Dim,T}) where {Dim,T} = new{Dim,T}(coords)
end

# convenience constructors
Point{Dim,T}(coords...) where {Dim,T} = Point{Dim,T}(coords)
Point{Dim,T}(coords) where {Dim,T} = Point{Dim,T}(Vec{Dim,T}(coords))
Point{Dim,T}(coords) where {Dim,T<:Integer} = Point{Dim,Float64}(coords)

Point(coords...) = Point(coords)
Point(coords) = Point(Vec(coords))
Point(coords::Vec{Dim,T}) where {Dim,T} = Point{Dim,T}(coords)

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

# broadcast behavior
Broadcast.broadcastable(p::Point) = Ref(p)

"""
    embeddim(point)

Return the number of dimensions of the space where the `point` is embedded.
"""
embeddim(::Type{Point{Dim,T}}) where {Dim,T} = Dim
embeddim(p::Point) = embeddim(typeof(p))

"""
    paramdim(point)

Return the number of parametric dimensions of the `point`.
"""
paramdim(::Type{Point{Dim,T}}) where {Dim,T} = 0
paramdim(p::Point) = paramdim(typeof(p))

"""
    coordtype(point)

Return the machine type of each coordinate used to describe the `point`.
"""
coordtype(::Type{Point{Dim,T}}) where {Dim,T} = T
coordtype(p::Point) = coordtype(typeof(p))

"""
    coordinates(A::Point)

Return the coordinates of the point with respect to the
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
    isapprox(A::Point, B::Point)

Tells whether or not the coordinates of points `A` and `B`
are approximately equal.
"""
Base.isapprox(A::Point{Dim,T}, B::Point{Dim,T}; atol=atol(T), kwargs...) where {Dim,T} =
  isapprox(A.coords, B.coords; atol, kwargs...)

"""
    ==(A::Point, B::Point)

Tells whether or not points `A` and `B` represent the same point.
"""
==(A::Point, B::Point) = A.coords == B.coords

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
    center(point)

Return the `point` itself.
"""
center(p::Point) = p

"""
    centroid(point)

Return the `point` itself.
"""
centroid(p::Point) = p

"""
    measure(point)

Return the measure of `point`, which is zero.
"""
measure(::Point{Dim,T}) where {Dim,T} = zero(T)

"""
    boundary(point)

Return the boundary of the `point`.
"""
boundary(::Point) = nothing

"""
    rand(P::Type{<:Point}, n=1)

Generates a random point of type `P`
"""
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Point{Dim,T}}) where {Dim,T} = Point(rand(rng, Vec{Dim,T}))

function Base.show(io::IO, point::Point)
  print(io, "Point$(Tuple(point.coords))")
end
