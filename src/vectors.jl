# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Vec(x₁, x₂, ..., xₙ)
    Vec((x₁, x₂, ..., xₙ))

A geometric vector in `Dim`-dimensional space with coordinates
in length units (default to meters) for linear algebra.

By default, integer coordinates are converted to float.

A vector can be obtained by subtracting two [`Point`](@ref) objects.

## Examples

```julia
A = Point(0.0, 0.0)
B = Point(1.0, 0.0)
v = B - A

# 2D vectors
Vec(1.0, 2.0) # add default units
Vec(1.0m, 2.0m) # double precision as expected
Vec(1f0km, 2f0km) # single precision as expected
Vec(1m, 2m) # integer is converted to float by design

# 3D vectors
Vec(1.0, 2.0, 3.0) # add default units
Vec(1.0m, 2.0m, 3.0m) # double precision as expected
Vec(1f0km, 2f0km, 3f0km) # single precision as expected
Vec(1m, 2m, 3m) # integer is converted to float by design
```

### Notes

A `Vec` is a subtype of `StaticVector` from StaticArrays.jl.
"""
struct Vec{Dim,ℒ<:Len} <: StaticVector{Dim,ℒ}
  coords::NTuple{Dim,ℒ}
  Vec{Dim,ℒ}(coords::NTuple{Dim}) where {Dim,ℒ<:Len} = new(coords)
end

Vec(coords::NTuple{Dim,ℒ}) where {Dim,ℒ<:Len} = Vec{Dim,float(ℒ)}(coords)
Vec(coords::NTuple{Dim,Len}) where {Dim} = Vec(promote(coords...))
Vec(coords::NTuple{Dim,Number}) where {Dim} = Vec(aslen.(coords))

Vec(coords::Number...) = Vec(coords)

# StaticVector interface
Base.Tuple(v::Vec) = getfield(v, :coords)
Base.getindex(v::Vec, i::Int) = getindex(getfield(v, :coords), i)
Base.promote_rule(::Type{Vec{Dim,ℒ₁}}, ::Type{Vec{Dim,ℒ₂}}) where {Dim,ℒ₁,ℒ₂} = Vec{Dim,promote_type(ℒ₁, ℒ₂)}
function StaticArrays.similar_type(::Type{<:Vec}, ::Type{T}, ::Size{S}) where {T,S}
  L = prod(S)
  N = length(S)
  isone(N) && T <: Len ? Vec{L,T} : SArray{Tuple{S...},T,N,L}
end

"""
    ∠(u, v)

Angle between vectors `u` and `v`.

Uses the two-argument form of `atan` returning
value in range [-π, π] in 2D and [0, π] in 3D.
See <https://en.wikipedia.org/wiki/Atan2>.

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```
"""
function ∠(u::Vec{2}, v::Vec{2}) # preserve sign
  θ = atan(u × v, u ⋅ v) * u"rad"
  θ == oftype(θ, -π) ? -θ : θ
end

∠(u::Vec{3}, v::Vec{3}) = atan(norm(u × v), u ⋅ v) * u"rad" # discard sign

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::Vec)
  if get(io, :compact, false)
    print(io, v.coords)
  else
    print(io, "Vec$(v.coords)")
  end
end

Base.show(io::IO, ::MIME"text/plain", v::Vec) = show(io, v)
