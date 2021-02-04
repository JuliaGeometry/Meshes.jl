# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Point{Dim,T}

A point in `Dim`-dimensional space with coordinates of type `T`.
The coordinates of the point provided upon construction are with
respect to the canonical Euclidean basis.

## Example

```julia
O = Point(0.0, 0.0) # origin of 2D Euclidean space
```

### Notes

- Type aliases are `Point1`, `Point2`, `Point3`, `Point1f`, `Point2f`, `Point3f`
"""
struct Point{Dim,T}
  coords::SVector{Dim,T}
  Point{Dim,T}(coords::AbstractVector) where {Dim,T} = new{Dim,T}(coords)
end

# convenience constructors
Point{Dim,T}(coords::NTuple{Dim,V}) where {Dim,T,V} = Point{Dim,T}(SVector(coords))
Point{Dim,T}(coords...) where {Dim,T} = Point{Dim,T}(coords)
Point(coords::NTuple{Dim,T}) where {Dim,T} = Point{Dim,T}(SVector(coords))
Point(coords...) = Point(coords)
Point(coords::AbstractVector{T}) where {T} =
  Point{length(coords),T}(SVector{length(coords)}(coords))

# catches mismatching tuple dimensions
Point{Dim1,T}(coords::NTuple{Dim2}) where {T,Dim1,Dim2} = 
  throw(DimensionMismatch("Can't construct a Point{$(Dim1),$(T)} with an input of length $(Dim2)"))

# coordinate type conversions
Base.convert(::Type{Point{Dim,T}}, coords) where {Dim,T} = Point{Dim,T}(coords)
Base.convert(::Type{Point{Dim,T}}, p::Point) where {Dim,T} = Point{Dim,T}(p.coords)
Base.convert(::Type{Point}, coords) = Point{length(coords),eltype(coords)}(coords)

# type aliases for convenience
const Point1 = Point{1,Float64}
const Point2  = Point{2,Float64}
const Point3  = Point{3,Float64}
const Point1f = Point{1,Float32}
const Point2f = Point{2,Float32}
const Point3f = Point{3,Float32}

"""
    embeddim(point)

Return the number of dimensions of the space where the `point` is embedded.
"""
embeddim(::Type{Point{Dim,T}}) where {Dim,T} = Dim
embeddim(p::Point) = embeddim(typeof(p))

"""
    coordtype(point)

Return the machine type of each coordinate used to describe the `point`.
"""
coordtype(::Type{Point{Dim,T}}) where {Dim,T}  = T
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
from point `A` to point `B`.
"""
-(A::Point, B::Point) = Vec(A.coords - B.coords)

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
Base.isapprox(A::Point, B::Point; kwargs...) = isapprox(A.coords, B.coords; kwargs...)

"""
    rand(P::Type{<:Point}, n=1)

Generates a random point of type `P`
"""
Random.rand(rng::Random.AbstractRNG,
            ::Random.SamplerType{Point{Dim,T}}) where {Dim,T} =
  Point(rand(rng, SVector{Dim,T}))

function Base.show(io::IO, point::Point)
  print(io, "Point$(Tuple(point.coords))")
end
