# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Vec(x₁, x₂, ..., xₙ)
    Vec((x₁, x₂, ..., xₙ))
    Vec([x₁, x₂, ..., xₙ])
    Vec{Dim,T}(x₁, x₂, ..., xₙ)
    Vec{Dim,T}((x₁, x₂, ..., xₙ))
    Vec{Dim,T}([x₁, x₂, ..., xₙ])

A vector in `Dim`-dimensional space with coordinates of type `T`.

By default, integer coordinates are converted to Float64.

A vector can be obtained by subtracting two [`Point`](@ref) objects.

## Examples

```julia
A = Point(0.0, 0.0)
B = Point(1.0, 0.0)
v = B - A

# 2D vectors
Vec(0.0, 1.0) # double precision as expected
Vec(0f0, 1f0) # single precision as expected
Vec(0, 0) # Integer is converted to Float64 by design
Vec2(0, 1) # explicitly ask for double precision
Vec2f(0, 1) # explicitly ask for single precision

# 3D vectors
Vec(1.0, 2.0, 3.0) # double precision as expected
Vec(1f0, 2f0, 3f0) # single precision as expected
Vec(1, 2, 3) # Integer is converted to Float64 by design
Vec3(1, 2, 3) # explicitly ask for double precision
Vec3f(1, 2, 3) # explicitly ask for single precision
```

### Notes

- A `Vec` is a subtype of `StaticVector` from StaticArrays.jl
- Type aliases are `Vec1`, `Vec2`, `Vec3`, `Vec1f`, `Vec2f`, `Vec3f`
"""
struct Vec{Dim,T} <: StaticVector{Dim,T}
  coords::SVector{Dim,T}
  Vec{Dim,T}(coords::SVector{Dim,T}) where {Dim,T} = new{Dim,T}(coords)
  Vec{Dim,T}(coords::SVector{Dim,T}) where {Dim,T<:Integer} = new{Dim,Float64}(coords)
end

# convenience constructors
Vec{Dim,T}(coords...) where {Dim,T} = Vec{Dim,T}(SVector{Dim,T}(coords))
Vec{Dim,T}(coords::Tuple) where {Dim,T} = Vec{Dim,T}(SVector{Dim,T}(coords))
Vec{Dim,T}(coords::AbstractVector) where {Dim,T} = Vec{Dim,T}(SVector{Dim,T}(coords))
Vec(coords::NTuple{Dim,T}) where {Dim,T} = Vec{Dim,T}(coords)
Vec(coords::SVector{Dim,T}) where {Dim,T} = Vec{Dim,T}(coords)
Vec(coords::AbstractVector{T}) where {T} = Vec{length(coords),T}(coords)
Vec(coords...) = Vec(SVector(coords))

# type aliases for convenience
const Vec1  = Vec{1,Float64}
const Vec2  = Vec{2,Float64}
const Vec3  = Vec{3,Float64}
const Vec1f = Vec{1,Float32}
const Vec2f = Vec{2,Float32}
const Vec3f = Vec{3,Float32}

# StaticVector interface
Base.Tuple(v::Vec) = Tuple(v.coords)
Base.getindex(v::Vec, i::Int) = getindex(v.coords, i)
function StaticArrays.similar_type(::Type{<:Vec}, ::Type{T}, ::Size{S}) where {T,S}
  L = prod(S)
  N = length(S)
  N == 1 ? Vec{L,T} : SArray{Tuple{S...}, T, N, L}
end
