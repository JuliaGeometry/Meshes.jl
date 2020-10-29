# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Point{Dim,T}

A point in `Dim`-dimensional space with coordinates of type `T`.
The coordinates of the point provided upon construction are with
respect to the canonical Euclidean basis. See [`vunit`](@ref).

## Example

```julia
O = Point(0.0, 0.0) # origin of 2D Euclidean space
```

### Notes

- Type aliases are `Point2`, `Point3`, `Point2f`, `Point3f`
"""
struct Point{Dim,T}
  coords::SVector{Dim,T}
end

# convenience constructors
Point(coords::NTuple{Dim,T}) where {Dim,T} = Point{Dim,T}(SVector(coords))
Point(coords::Vararg{T,Dim}) where {Dim,T} = Point{Dim,T}(SVector(coords))

# coordinate type conversions
Point{Dim,T}(coords::NTuple{Dim,V}) where {Dim,T,V} = Point(T.(coords))
Point{Dim,T}(coords::Vararg{V,Dim}) where {Dim,T,V} = Point(T.(coords))
Base.convert(::Type{Point{Dim,T}}, coords) where {Dim,T} = Point{Dim,T}(coords)
Base.convert(::Type{Point{Dim,T}}, p::Point) where {Dim,T} = Point{Dim,T}(p.coords)
Base.convert(::Type{Point}, coords) = Point{length(coords),eltype(coords)}(coords)

# type aliases for convenience
const Point2  = Point{2,Float64}
const Point3  = Point{3,Float64}
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

"""
    -(A::Point, v::Vec)
    -(v::Vec, A::Point)

Return the point at the end of the vector `-v` placed
at a reference (or start) point `A`.
"""
-(A::Point, v::Vec) = Point(A.coords - v)
-(v::Vec, A::Point) = A - v

# TODO: implement rand properly with RNG, etc.
Base.rand(::Type{Point{Dim,T}}) where {Dim,T} = Point(rand(SVector{Dim,T}))
Base.rand(::Type{Point{Dim,T}}, n::Integer) where {Dim,T} = Point.(rand(SVector{Dim,T}, n))

# -----------
# IO methods
# -----------
function Base.show(io::IO, point::Point)
  print(io, "Point$(Tuple(point.coords))")
end
