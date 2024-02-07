# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cartesian(x₁, x₂, ..., xₙ)

N-dimensional Cartesian coordinates `x₁, x₂, ..., xₙ` in length units (default to meter).

## Examples

```julia
Cartesian(1, 1) # add default units
Cartesian(1u"m", 1u"m") # integers are converted converted to floats
Cartesian(1.0u"km", 1.0u"km", 1.0u"km")
```

## References

* [Cartesian coordinate system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
"""
struct Cartesian{N,L<:Len} <: Coordinates{N}
  coords::NTuple{N,L}
  Cartesian{N,L}(coords) where {N,L} = new{N,float(L)}(coords)
end

Cartesian(coords::Vararg{L,N}) where {N,L<:Len} = Cartesian{N,L}(coords)
Cartesian(coords::Len...) = Cartesian(promote(coords...)...)
Cartesian(coords::Number...) = Cartesian(addunit.(coords, u"m")...)

_fields(coords::Cartesian) = coords.coords
function _fnames(::Cartesian{N}) where {N}
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

Polar coordinates with radius `ρ ∈ [0,∞)` in length units (default to meter)
and angle `ϕ ∈ [0,2π)` in angular units (default to radian).

## Examples

```julia
Polar(1, π/4) # add default units
Polar(1u"m", (π/4)u"rad") # integers are converted converted to floats
Polar(1.0u"m", 45u"°") # degrees are converted to radians
Polar(1.0u"km", (π/4)u"rad")
```

## References

* [Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
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

Cylindrical coordinates with radius `ρ ∈ [0,∞)` in length units (default to meter), 
angle `ϕ ∈ [0,2π)` in angular units (default to radian) 
and height `z ∈ [0,∞)` in length units (default to meter).

## Examples

```julia
Cylindrical(1, π/4, 1) # add default units
Cylindrical(1u"m", (π/4)u"rad", 1u"m") # integers are converted converted to floats
Cylindrical(1.0u"m", 45u"°", 1.0u"m") # degrees are converted to radians
Cylindrical(1.0u"km", (π/4)u"rad", 1.0u"km")
```

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
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

Spherical coordinates with radius `r ∈ [0,∞)` in length units (default to meter), 
polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)` in angular units (default to radian).

## Examples

```julia
Spherical(1, π/4, π/4) # add default units
Spherical(1u"m", (π/4)u"rad", (π/4)u"rad") # integers are converted converted to floats
Spherical(1.0u"m", 45u"°", 45u"°") # degrees are converted to radians
Spherical(1.0u"km", (π/4)u"rad", (π/4)u"rad")
```

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
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

# ------------
# CONVERSIONS
# ------------

# Cartesian <> Polar
Base.convert(::Type{Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
function Base.convert(::Type{Polar}, (; coords)::Cartesian{2})
  x, y = coords
  Polar(sqrt(x^2 + y^2), atanpos(y, x) * u"rad")
end

# Cartesian <> Cylindrical
Base.convert(::Type{Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
function Base.convert(::Type{Cylindrical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Cylindrical(sqrt(x^2 + y^2), atanpos(y, x) * u"rad", z)
end

# Cartesian <> Spherical
Base.convert(::Type{Cartesian}, (; r, θ, ϕ)::Spherical) =
  Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
function Base.convert(::Type{Spherical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z) * u"rad", atanpos(y, x) * u"rad")
end
