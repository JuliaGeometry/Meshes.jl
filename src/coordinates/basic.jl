# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cartesian(x₁, x₂, ..., xₙ)

N-dimensional Cartesian coordinates `x₁, x₂, ..., xₙ`.

## References

* [Cartesian coordinate system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Cartesian{N,T<:Len} <: Coordinates{N}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::Vararg{T,N}) where {N,T<:Len} = Cartesian{N,T}(coords)
Cartesian(coords::Len...) = Cartesian(promote(coords...)...)
Cartesian(coords::Number...) = Cartesian((coords .* u"m")...)

Base.isapprox(c₁::Cartesian{N}, c₂::Cartesian{N}; kwargs...) where {N} =
  all(isapprox(x₁, x₂; kwargs...) for (x₁, x₂) in zip(c₁.coords, c₂.coords))

function Base.show(io::IO, (; coords)::Cartesian{N}) where {N}
  print(io, "Cartesian(")
  fnames = _cartfields(N)
  printfields(io, coords, fnames, compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", (; coords)::Cartesian{N}) where {N}
  print(io, "Cartesian coordinates")
  fnames = _cartfields(N)
  printfields(io, coords, fnames)
end

function _cartfields(N)
  if N == 1
    ("x",)
  elseif N == 2
    ("x", "y")
  elseif N == 3
    ("x", "y", "z")
  else
    ntuple(i -> "x$i", N)
  end
end

"""
    Polar(ρ, ϕ)

Polar coordinates with radius `ρ ∈ [0,∞)` and angle `ϕ ∈ [0,2π)`.

## References

* [Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Polar{T<:Len,A<:Rad} <: Coordinates{2}
  ρ::T
  ϕ::A
  Polar{T,A}(ρ, ϕ) where {T,A} = new{float(T),float(A)}(ρ, ϕ)
end

Polar(ρ::T, ϕ::A) where {T<:Len,A<:Rad} = Polar{T,A}(ρ, ϕ)
Polar(ρ::Len, ϕ::Deg) = Polar(ρ, deg2rad(ϕ))
Polar(ρ::Number, ϕ::Number) = Polar(ρ * u"m", ϕ * u"rad")

"""
    Cylindrical(ρ, ϕ, z)

Cylindrical coordinates with radius `ρ ∈ [0,∞)`, angle `ϕ ∈ [0,2π)` and height `z ∈ [0,∞)`.

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Cylindrical{T<:Len,A<:Rad} <: Coordinates{3}
  ρ::T
  ϕ::A
  z::T
  Cylindrical{T,A}(ρ, ϕ, z) where {T,A} = new{float(T),float(A)}(ρ, ϕ, z)
end

Cylindrical(ρ::T, ϕ::A, z::T) where {T<:Len,A<:Rad} = Cylindrical{T,A}(ρ, ϕ, z)
function Cylindrical(ρ::Len, ϕ::Rad, z::Len)
  nρ, nz = promote(ρ, z)
  Cylindrical(nρ, ϕ, nz)
end
Cylindrical(ρ::Len, ϕ::Deg, z::Len) = Cylindrical(ρ, deg2rad(ϕ), z)
Cylindrical(ρ::Number, ϕ::Number, z::Number) = Cylindrical(ρ * u"m", ϕ * u"rad", z * u"m")

"""
    Spherical(r, θ, ϕ)

Spherical coordinates with radius `r ∈ [0,∞)`, polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)`.

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Spherical{T<:Len,A<:Rad} <: Coordinates{3}
  r::T
  θ::A
  ϕ::A
  Spherical{T,A}(r, θ, ϕ) where {T,A} = new{float(T),float(A)}(r, θ, ϕ)
end

Spherical(r::T, θ::A, ϕ::A) where {T<:Len,A<:Rad} = Spherical{T,A}(r, θ, ϕ)
Spherical(r::Len, θ::Rad, ϕ::Rad) = Spherical(r, promote(θ, ϕ)...)
Spherical(r::Len, θ::Deg, ϕ::Deg) = Spherical(r, deg2rad(θ), deg2rad(ϕ))
Spherical(r::Number, θ::Number, ϕ::Number) = Spherical(r * u"m", θ * u"rad", ϕ * u"rad")
