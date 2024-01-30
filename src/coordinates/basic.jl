# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Cartesian{N,T} <: Coordinates{N,T}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::NTuple{N,T}) where {N,T} = Coordinates{N,T}(coords)
Cartesian(coords::Tuple) = Coordinates(promote(coords...))
Cartesian(coords...) = Cartesian(coords)

struct Polar{T} <: Coordinates{2,T}
  ρ::T
  ϕ::T
  Polar{T}(ρ, ϕ) where {T} = new{float(T)}(ρ, ϕ)
end

Polar(ρ::T, ϕ::T) where {T} = Polar{T}(ρ, ϕ)
Polar(ρ, ϕ) = Polar(promote(ρ, ϕ)...)

struct Cylindrical{T} <: Coordinates{3,T}
  ρ::T
  ϕ::T
  z::T
  Cylindrical{T}(ρ, ϕ, z) where {T} = new{float(T)}(ρ, ϕ, z)
end

Cylindrical(ρ::T, ϕ::T, z::T) where {T} = Cylindrical{T}(ρ, ϕ, z)
Cylindrical(ρ, ϕ, z) = Cylindrical(promote(ρ, ϕ, z)...)

struct Spherical{T} <: Coordinates{3,T}
  r::T
  θ::T
  ϕ::T
  Spherical{T}(r, θ, ϕ) where {T} = new{float(T)}(r, θ, ϕ)
end

Spherical(r::T, θ::T, ϕ::T) where {T} = Spherical{T}(r, θ, ϕ)
Spherical(r, θ, ϕ) = Spherical(promote(r, θ, ϕ)...)
