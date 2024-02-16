# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cartesian(x₁, x₂, ..., xₙ)

N-dimensional Cartesian coordinates `x₁, x₂, ..., xₙ` in length units (default to meter).
The first 3 coordinates can be accessed with the properties `x`, `y` and `z`, respectively.

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
struct Cartesian{N,L<:Len} <: CRS{NoDatum}
  coords::NTuple{N,L}
  Cartesian(coords::NTuple{N,L}) where {N,L<:Len} = new{N,float(L)}(coords)
end

Cartesian(coords::L...) where {L<:Len} = Cartesian(coords)
Cartesian(coords::Len...) = Cartesian(promote(coords...))
Cartesian(coords::Number...) = Cartesian(addunit.(coords, u"m")...)

Base.propertynames(::Cartesian) = (:x, :y, :z)

function Base.getproperty(coords::Cartesian, name::Symbol)
  tup = _coords(coords)
  if name === :x
    tup[1]
  elseif name === :y
    tup[2]
  elseif name === :z
    tup[3]
  else
    error("invalid property, use `x`, `y` or `z`")
  end
end

_coords(coords::Cartesian) = getfield(coords, :coords)

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
struct Polar{L<:Len,R<:Rad} <: CRS{NoDatum}
  ρ::L
  ϕ::R
  Polar(ρ::L, ϕ::R) where {L<:Len,R<:Rad} = new{float(L),float(R)}(ρ, ϕ)
end

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
struct Cylindrical{L<:Len,R<:Rad} <: CRS{NoDatum}
  ρ::L
  ϕ::R
  z::L
  Cylindrical(ρ::L, ϕ::R, z::L) where {L<:Len,R<:Rad} = new{float(L),float(R)}(ρ, ϕ, z)
end

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
struct Spherical{L<:Len,R<:Rad} <: CRS{NoDatum}
  r::L
  θ::R
  ϕ::R
  Spherical(r::L, θ::R, ϕ::R) where {L<:Len,R<:Rad} = new{float(L),float(R)}(r, θ, ϕ)
end

Spherical(r::Len, θ::Rad, ϕ::Rad) = Spherical(r, promote(θ, ϕ)...)
Spherical(r::Len, θ::Deg, ϕ::Deg) = Spherical(r, deg2rad(θ), deg2rad(ϕ))
Spherical(r::Number, θ::Number, ϕ::Number) = Spherical(addunit(r, u"m"), addunit(θ, u"rad"), addunit(ϕ, u"rad"))

# ------------
# CONVERSIONS
# ------------

# Cartesian <> Polar
Base.convert(::Type{Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
Base.convert(::Type{Polar}, (; x, y)::Cartesian{2}) = Polar(sqrt(x^2 + y^2), atanpos(y, x) * u"rad")

# Cartesian <> Cylindrical
Base.convert(::Type{Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
Base.convert(::Type{Cylindrical}, (; x, y, z)::Cartesian{3}) = Cylindrical(sqrt(x^2 + y^2), atanpos(y, x) * u"rad", z)

# Cartesian <> Spherical
Base.convert(::Type{Cartesian}, (; r, θ, ϕ)::Spherical) =
  Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
Base.convert(::Type{Spherical}, (; x, y, z)::Cartesian{3}) =
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z) * u"rad", atanpos(y, x) * u"rad")
