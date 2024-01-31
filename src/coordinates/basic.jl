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
struct Cartesian{N,T} <: Coordinates{N}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::NTuple{N,T}) where {N,T} = Cartesian{N,T}(coords)
Cartesian(coords::Tuple) = Cartesian(promote(coords...))
Cartesian(coords...) = Cartesian(coords)

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
struct Polar{T,A} <: Coordinates{2}
  ρ::T
  ϕ::A
  Polar{T,A}(ρ, ϕ) where {T,A} = new{float(T),float(A)}(ρ, ϕ)
end

Polar(ρ::T, ϕ::A) where {T,A} = Polar{T,A}(ρ, ϕ)

"""
    Cylindrical(ρ, ϕ, z)

Cylindrical coordinates with radius `ρ ∈ [0,∞)`, angle `ϕ ∈ [0,2π)` and height `z ∈ [0,∞)`.

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Cylindrical{T,A} <: Coordinates{3}
  ρ::T
  ϕ::A
  z::T
  Cylindrical{T,A}(ρ, ϕ, z) where {T,A} = new{float(T),float(A)}(ρ, ϕ, z)
end

Cylindrical(ρ::T, ϕ::A, z::T) where {T,A} = Cylindrical{T,A}(ρ, ϕ, z)
function Cylindrical(ρ, ϕ, z)
  nρ, nz = promote(ρ, z)
  Cylindrical(nρ, ϕ, nz)
end

"""
    Spherical(r, θ, ϕ)

Spherical coordinates with radius `r ∈ [0,∞)`, polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)`.

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Spherical{T,A} <: Coordinates{3}
  r::T
  θ::A
  ϕ::A
  Spherical{T,A}(r, θ, ϕ) where {T,A} = new{float(T),float(A)}(r, θ, ϕ)
end

Spherical(r::T, θ::A, ϕ::A) where {T,A} = Spherical{T,A}(r, θ, ϕ)
Spherical(r, θ, ϕ) = Spherical(r, promote(θ, ϕ)...)
