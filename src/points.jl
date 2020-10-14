abstract type AbstractPoint{N,T} end

"""
    Point{N,T}

A point in `N`-dimensional space with coordinates of type `T`.
The coordinates of the point provided upon construction are with
respect to the canonical Euclidean basis. See [`vunit`](@ref).

## Example

```julia
O = Point(0.0, 0.0) # origin of 2D Euclidean space
```

### Notes

- Type aliases are `Point2`, `Point3`, `Point2f`, `Point3f`
"""
struct Point{N,T} <: AbstractPoint{N,T}
    coords::SVector{N,T}
end

# convenience constructors
Point(coords::NTuple{N,T}) where {N,T} = Point{N,T}(SVector(coords))
Point(coords::Vararg{T,N}) where {N,T} = Point{N,T}(SVector(coords))

# coordinate type conversions
Point{N,T}(coords::NTuple{N,V}) where {N,T,V} = Point(T.(coords))
Point{N,T}(coords::Vararg{V,N}) where {N,T,V} = Point(T.(coords))
Base.convert(::Type{Point{N,T}}, coords) where {N,T} = Point{N,T}(coords)
Base.convert(::Type{Point{N,T}}, p::Point) where {N,T} = Point{N,T}(p.coords)

# type aliases for convenience
const Point2  = Point{2,Float64}
const Point3  = Point{3,Float64}
const Point2f = Point{2,Float32}
const Point3f = Point{3,Float32}

"""
    coordinates(A::Point)

Return the coordinates of the point with respect to the
canonical Euclidean basis. See [`vunit`](@ref).
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

# TODO: implement rand properly with RNG, etc.
Base.rand(::Type{Point{N,T}}) where {N,T} = Point(rand(SVector{N,T}))
Base.rand(::Type{Point{N,T}}, n::Integer) where {N,T} = Point.(rand(SVector{N,T}, n))
