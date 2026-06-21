# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Vec(x, [y, z])
    Vec((x, [y, z]))

A vector with `Cartesian` coordinates in length units (default to meters).
Vectors are used in linear algebra and can be obtained by subtracting two
[`Point`](@ref)s in any coordinate reference system.

See [`Point`](@ref) for more details on preprocessing of coordinates and units.

## Examples

```julia
# vectors in 2D Euclidean space
Vec(1.0, 2.0)
Vec(1km, 2km)

# vectors in 3D Euclidean space
Vec(1.0, 2.0, 3.0)
Vec(1km, 2km, 3km)

# vector as the difference of two points
Point(1.0, 1.0) - Point(0.0, 0.0)
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
