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
struct Cartesian{N,T} <: Coordinates{N,T}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::NTuple{N,T}) where {N,T} = Cartesian{N,T}(coords)
Cartesian(coords::Tuple) = Cartesian(promote(coords...))
Cartesian(coords...) = Cartesian(coords)

"""
    Polar(ρ, ϕ)

Polar coordinates with radius `ρ ∈ [0,∞)` and angle `ϕ ∈ [0,2π)`.

## References

* [Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Polar{T} <: Coordinates{2,T}
  ρ::T
  ϕ::T
  Polar{T}(ρ, ϕ) where {T} = new{float(T)}(ρ, ϕ)
end

Polar(ρ::T, ϕ::T) where {T} = Polar{T}(ρ, ϕ)
Polar(ρ, ϕ) = Polar(promote(ρ, ϕ)...)

"""
    Cylindrical(ρ, ϕ, z)

Cylindrical coordinates with radius `ρ ∈ [0,∞)`, angle `ϕ ∈ [0,2π)` and height `z ∈ [0,∞)`.

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Cylindrical{T} <: Coordinates{3,T}
  ρ::T
  ϕ::T
  z::T
  Cylindrical{T}(ρ, ϕ, z) where {T} = new{float(T)}(ρ, ϕ, z)
end

Cylindrical(ρ::T, ϕ::T, z::T) where {T} = Cylindrical{T}(ρ, ϕ, z)
Cylindrical(ρ, ϕ, z) = Cylindrical(promote(ρ, ϕ, z)...)

"""
    Spherical(r, θ, ϕ)

Spherical coordinates with radius `r ∈ [0,∞)`, polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)`.

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Spherical{T} <: Coordinates{3,T}
  r::T
  θ::T
  ϕ::T
  Spherical{T}(r, θ, ϕ) where {T} = new{float(T)}(r, θ, ϕ)
end

Spherical(r::T, θ::T, ϕ::T) where {T} = Spherical{T}(r, θ, ϕ)
Spherical(r, θ, ϕ) = Spherical(promote(r, θ, ϕ)...)
