# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cartesian(x₁, x₂, ..., xₙ)
    Cartesian{Datum}(x₁, x₂, ..., xₙ)

N-dimensional Cartesian coordinates `x₁, x₂, ..., xₙ` in length units (default to meter)
with a given `Datum` (default to `NoDatum`).
The first 3 coordinates can be accessed with the properties `x`, `y` and `z`, respectively.

## Examples

```julia
Cartesian(1, 1) # add default units
Cartesian(1u"m", 1u"m") # integers are converted converted to floats
Cartesian(1.0u"km", 1.0u"km", 1.0u"km")
Cartesian{WGS84}(1.0u"m", 1.0u"m")
```

## References

* [Cartesian coordinate system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
"""
struct Cartesian{Datum,N,L<:Len} <: CRS{Datum}
  coords::NTuple{N,L}
  Cartesian{Datum}(coords::NTuple{N,L}) where {Datum,N,L<:Len} = new{Datum,N,float(L)}(coords)
end

Cartesian{Datum}(coords::L...) where {Datum,L<:Len} = Cartesian{Datum}(coords)
Cartesian{Datum}(coords::Len...) where {Datum} = Cartesian{Datum}(promote(coords...))
Cartesian{Datum}(coords::Number...) where {Datum} = Cartesian{Datum}(addunit.(coords, u"m")...)

Cartesian(args...) = Cartesian{NoDatum}(args...)

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

function Base.isapprox(coords₁::C, coords₂::C; kwargs...) where {C<:Cartesian}
  tup₁ = _coords(coords₁)
  tup₂ = _coords(coords₂)
  all(isapprox(c₁, c₂; kwargs...) for (c₁, c₂) in zip(tup₁, tup₂))
end

function Base.show(io::IO, coords::Cartesian)
  Datum = datum(coords)
  print(io, "Cartesian{$Datum}(")
  printfields(io, _coords(coords), _fnames(coords), compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", coords::Cartesian)
  Datum = datum(coords)
  print(io, "Cartesian{$Datum} coordinates")
  printfields(io, _coords(coords), _fnames(coords))
end

_coords(coords::Cartesian) = getfield(coords, :coords)

function _fnames(::Cartesian{<:Any,N}) where {N}
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
    Polar{Datum}(ρ, ϕ)

Polar coordinates with radius `ρ ∈ [0,∞)` in length units (default to meter),
angle `ϕ ∈ [0,2π)` in angular units (default to radian)
and a given `Datum` (default to `NoDatum`).

## Examples

```julia
Polar(1, π/4) # add default units
Polar(1u"m", (π/4)u"rad") # integers are converted converted to floats
Polar(1.0u"m", 45u"°") # degrees are converted to radians
Polar(1.0u"km", (π/4)u"rad")
Polar{WGS84}(1.0u"m", (π/4)u"rad")
```

## References

* [Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
"""
struct Polar{Datum,L<:Len,R<:Rad} <: CRS{Datum}
  ρ::L
  ϕ::R
  Polar{Datum}(ρ::L, ϕ::R) where {Datum,L<:Len,R<:Rad} = new{Datum,float(L),float(R)}(ρ, ϕ)
end

Polar{Datum}(ρ::Len, ϕ::Deg) where {Datum} = Polar{Datum}(ρ, deg2rad(ϕ))
Polar{Datum}(ρ::Number, ϕ::Number) where {Datum} = Polar{Datum}(addunit(ρ, u"m"), addunit(ϕ, u"rad"))

Polar(args...) = Polar{NoDatum}(args...)

"""
    Cylindrical(ρ, ϕ, z)
    Cylindrical{Datum}(ρ, ϕ, z)

Cylindrical coordinates with radius `ρ ∈ [0,∞)` in length units (default to meter), 
angle `ϕ ∈ [0,2π)` in angular units (default to radian),
height `z ∈ [0,∞)` in length units (default to meter)
and a given `Datum` (default to `NoDatum`).

## Examples

```julia
Cylindrical(1, π/4, 1) # add default units
Cylindrical(1u"m", (π/4)u"rad", 1u"m") # integers are converted converted to floats
Cylindrical(1.0u"m", 45u"°", 1.0u"m") # degrees are converted to radians
Cylindrical(1.0u"km", (π/4)u"rad", 1.0u"km")
Cylindrical{WGS84}(1.0u"m", (π/4)u"rad", 1.0u"m")
```

## References

* [Cylindrical coordinate system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
"""
struct Cylindrical{Datum,L<:Len,R<:Rad} <: CRS{Datum}
  ρ::L
  ϕ::R
  z::L
  Cylindrical{Datum}(ρ::L, ϕ::R, z::L) where {Datum,L<:Len,R<:Rad} = new{Datum,float(L),float(R)}(ρ, ϕ, z)
end

function Cylindrical{Datum}(ρ::Len, ϕ::Rad, z::Len) where {Datum}
  nρ, nz = promote(ρ, z)
  Cylindrical{Datum}(nρ, ϕ, nz)
end
Cylindrical{Datum}(ρ::Len, ϕ::Deg, z::Len) where {Datum} = Cylindrical{Datum}(ρ, deg2rad(ϕ), z)
Cylindrical{Datum}(ρ::Number, ϕ::Number, z::Number) where {Datum} =
  Cylindrical{Datum}(addunit(ρ, u"m"), addunit(ϕ, u"rad"), addunit(z, u"m"))

Cylindrical(args...) = Cylindrical{NoDatum}(args...)

"""
    Spherical(r, θ, ϕ)
    Spherical{Datum}(r, θ, ϕ)

Spherical coordinates with radius `r ∈ [0,∞)` in length units (default to meter), 
polar angle `θ ∈ [0,π]` and azimuth angle `ϕ ∈ [0,2π)` in angular units (default to radian)
and a given `Datum` (default to `NoDatum`).

## Examples

```julia
Spherical(1, π/4, π/4) # add default units
Spherical(1u"m", (π/4)u"rad", (π/4)u"rad") # integers are converted converted to floats
Spherical(1.0u"m", 45u"°", 45u"°") # degrees are converted to radians
Spherical(1.0u"km", (π/4)u"rad", (π/4)u"rad")
Spherical{WGS84}(1.0u"m", (π/4)u"rad", (π/4)u"rad")
```

## References

* [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ISO 80000-2:2019](https://www.iso.org/standard/64973.html)
* [ISO 80000-3:2019](https://www.iso.org/standard/64974.html)
"""
struct Spherical{Datum,L<:Len,R<:Rad} <: CRS{Datum}
  r::L
  θ::R
  ϕ::R
  Spherical{Datum}(r::L, θ::R, ϕ::R) where {Datum,L<:Len,R<:Rad} = new{Datum,float(L),float(R)}(r, θ, ϕ)
end

Spherical{Datum}(r::Len, θ::Rad, ϕ::Rad) where {Datum} = Spherical{Datum}(r, promote(θ, ϕ)...)
Spherical{Datum}(r::Len, θ::Deg, ϕ::Deg) where {Datum} = Spherical{Datum}(r, deg2rad(θ), deg2rad(ϕ))
Spherical{Datum}(r::Number, θ::Number, ϕ::Number) where {Datum} =
  Spherical{Datum}(addunit(r, u"m"), addunit(θ, u"rad"), addunit(ϕ, u"rad"))

Spherical(args...) = Spherical{NoDatum}(args...)

# ------------
# CONVERSIONS
# ------------

# Cartesian <> Polar
Base.convert(::Type{Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
Base.convert(::Type{Polar}, (; x, y)::Cartesian{<:Any,2}) = Polar(sqrt(x^2 + y^2), atanpos(y, x) * u"rad")

# Cartesian <> Cylindrical
Base.convert(::Type{Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
Base.convert(::Type{Cylindrical}, (; x, y, z)::Cartesian{<:Any,3}) =
  Cylindrical(sqrt(x^2 + y^2), atanpos(y, x) * u"rad", z)

# Cartesian <> Spherical
Base.convert(::Type{Cartesian}, (; r, θ, ϕ)::Spherical) =
  Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
Base.convert(::Type{Spherical}, (; x, y, z)::Cartesian{<:Any,3}) =
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z) * u"rad", atanpos(y, x) * u"rad")
