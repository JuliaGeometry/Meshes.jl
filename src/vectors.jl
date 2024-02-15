# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Continuous = Union{AbstractFloat,Quantity{<:AbstractFloat}}

"""
    Vec(x₁, x₂, ..., xₙ)
    Vec((x₁, x₂, ..., xₙ))
    Vec{Dim,T}(x₁, x₂, ..., xₙ)
    Vec{Dim,T}((x₁, x₂, ..., xₙ))

A geometric vector in `Dim`-dimensional space with coordinates of type `T` for linear algebra.

By default, integer coordinates are converted to float.

A vector can be obtained by subtracting two [`Point`](@ref) objects.

## Examples

```julia
A = Point(0.0, 0.0)
B = Point(1.0, 0.0)
v = B - A

# 2D vectors
Vec(0.0, 1.0) # double precision as expected
Vec(0f0, 1f0) # single precision as expected
Vec(0, 0) # integer is converted to float by design
Vec2(0, 1) # explicitly ask for double precision
Vec2f(0, 1) # explicitly ask for single precision

# 3D vectors
Vec(1.0, 2.0, 3.0) # double precision as expected
Vec(1f0, 2f0, 3f0) # single precision as expected
Vec(1, 2, 3) # integer is converted to float by design
Vec3(1, 2, 3) # explicitly ask for double precision
Vec3f(1, 2, 3) # explicitly ask for single precision
```

### Notes

- A `Vec` is a subtype of `StaticVector` from StaticArrays.jl
- Type aliases are `Vec1`, `Vec2`, `Vec3`, `Vec1f`, `Vec2f`, `Vec3f`
"""
struct Vec{Dim,T<:Continuous} <: StaticVector{Dim,T}
  coords::NTuple{Dim,T}
  Vec{Dim,T}(coords::NTuple{Dim,T}) where {Dim,T<:Continuous} = new(coords)
end

# convenience constructors
Vec{Dim,T}(coords::Number...) where {Dim,T<:Continuous} = Vec{Dim,T}(coords)
function Vec{Dim,T}(coords::Union{Tuple,AbstractVector}) where {Dim,T<:Continuous}
  if Dim ≠ length(coords)
    throw(DimensionMismatch("the number of coordinates must be equal to the number of dimensions"))
  end
  Vec{Dim,T}(NTuple{Dim,T}(coords))
end

Vec(coords::Number...) = Vec(coords)
Vec(coords::NTuple{Dim,Number}) where {Dim} = Vec(promote(coords...))
Vec(coords::NTuple{Dim,T}) where {Dim,T<:Number} = Vec(float.(coords))
Vec(coords::NTuple{Dim,T}) where {Dim,T<:Continuous} = Vec{Dim,T}(coords)

# StaticVector constructors
Vec(coords::StaticVector{Dim,T}) where {Dim,T} = Vec{Dim,T}(coords)
Vec{Dim,T}(coords::StaticVector) where {Dim,T} = Vec{Dim,T}(Tuple(coords))

# type aliases for convenience
const Vec1 = Vec{1,Float64}
const Vec2 = Vec{2,Float64}
const Vec3 = Vec{3,Float64}
const Vec1f = Vec{1,Float32}
const Vec2f = Vec{2,Float32}
const Vec3f = Vec{3,Float32}

# StaticVector interface
Base.Tuple(v::Vec) = getfield(v, :coords)
Base.getindex(v::Vec, i::Int) = getindex(getfield(v, :coords), i)
function StaticArrays.similar_type(::Type{<:Vec}, ::Type{T}, ::Size{S}) where {T,S}
  L = prod(S)
  N = length(S)
  isone(N) && T <: Continuous ? Vec{L,T} : SArray{Tuple{S...},T,N,L}
end

"""
    ∠(u, v)

Angle between vectors `u` and `v`.

Uses the two-argument form of `atan` returning
value in range [-π, π] in 2D and [0, π] in 3D.
See https://en.wikipedia.org/wiki/Atan2.

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```
"""
function ∠(u::Vec{2}, v::Vec{2}) # preserve sign
  θ = atan(u × v, u ⋅ v)
  T = typeof(θ)
  θ == -T(π) ? -θ : θ
end

∠(u::Vec{3}, v::Vec{3}) = atan(norm(u × v), u ⋅ v) # discard sign

function Base.show(io::IO, v::Vec)
  if get(io, :compact, false)
    print(io, v.coords)
  else
    print(io, "Vec$(v.coords)")
  end
end

Base.show(io::IO, ::MIME"text/plain", v::Vec) = show(io, v)
