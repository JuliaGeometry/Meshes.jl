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
struct Cartesian{N,L<:Len} <: Coordinates{N}
  coords::NTuple{N,L}
  Cartesian{N,L}(coords) where {N,L} = new{N,float(L)}(coords)
end

Cartesian(coords::Vararg{L,N}) where {N,L<:Len} = Cartesian{N,L}(coords)
Cartesian(coords::Len...) = Cartesian(promote(coords...)...)
Cartesian(coords::Number...) = Cartesian(addunit.(coords, u"m")...)

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
struct Polar{L<:Len,R<:Rad} <: Coordinates{2}
  ρ::L
  ϕ::R
  Polar{L,R}(ρ, ϕ) where {L,R} = new{float(L),float(R)}(ρ, ϕ)
end

Polar(ρ::L, ϕ::R) where {L<:Len,R<:Rad} = Polar{L,R}(ρ, ϕ)
Polar(ρ::Len, ϕ::Deg) = Polar(ρ, deg2rad(ϕ))
Polar(ρ::Number, ϕ::Number) = Polar(addunit(ρ, u"m"), addunit(ϕ, u"rad"))

"""
    Cylindrical(ρ, ϕ, z)

Cylindrical coordinates with radius `ρ ∈ [0,∞)`, angle `ϕ ∈ [0,2π)` and height `z ∈ [0,∞)`.

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Cylindrical{L<:Len,R<:Rad} <: Coordinates{3}
  ρ::L
  ϕ::R
  z::L
  Cylindrical{L,R}(ρ, ϕ, z) where {L,R} = new{float(L),float(R)}(ρ, ϕ, z)
end

Cylindrical(ρ::L, ϕ::R, z::L) where {L<:Len,R<:Rad} = Cylindrical{L,R}(ρ, ϕ, z)
function Cylindrical(ρ::Len, ϕ::Rad, z::Len)
  nρ, nz = promote(ρ, z)
  Cylindrical(nρ, ϕ, nz)
end
Cylindrical(ρ::Len, ϕ::Deg, z::Len) = Cylindrical(ρ, deg2rad(ϕ), z)
Cylindrical(ρ::Number, ϕ::Number, z::Number) = Cylindrical(addunit(ρ, u"m"), addunit(ϕ, u"rad"), addunit(z, u"m"))

"""
    Spherical(r, θ, ϕ)

Spherical coordinates with radius `r ∈ [0,∞)`, polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)`.

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 31-11](https://en.wikipedia.org/wiki/ISO_31-11)
"""
struct Spherical{L<:Len,R<:Rad} <: Coordinates{3}
  r::L
  θ::R
  ϕ::R
  Spherical{L,R}(r, θ, ϕ) where {L,R} = new{float(L),float(R)}(r, θ, ϕ)
end

Spherical(r::L, θ::R, ϕ::R) where {L<:Len,R<:Rad} = Spherical{L,R}(r, θ, ϕ)
Spherical(r::Len, θ::Rad, ϕ::Rad) = Spherical(r, promote(θ, ϕ)...)
Spherical(r::Len, θ::Deg, ϕ::Deg) = Spherical(r, deg2rad(θ), deg2rad(ϕ))
Spherical(r::Number, θ::Number, ϕ::Number) = Spherical(addunit(r, u"m"), addunit(θ, u"rad"), addunit(ϕ, u"rad"))
