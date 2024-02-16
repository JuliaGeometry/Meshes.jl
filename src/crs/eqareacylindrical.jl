# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EqualAreaCylindrical{Datum,latₜₛ}

Equal Area Cylindrical CRS with a given `Datum` and latitude of true scale `latₜₛ` in degrees.
"""
struct EqualAreaCylindrical{Datum,latₜₛ,M<:Met} <: CRS{Datum}
  x::M
  y::M
  EqualAreaCylindrical{Datum,latₜₛ}(x::M, y::M) where {Datum,latₜₛ,M<:Met} = new{Datum,latₜₛ,float(M)}(x, y)
end

EqualAreaCylindrical{Datum,latₜₛ}(x::Met, y::Met) where {Datum,latₜₛ} =
  EqualAreaCylindrical{Datum,latₜₛ}(promote(x, y)...)
EqualAreaCylindrical{Datum,latₜₛ}(x::Len, y::Len) where {Datum,latₜₛ} =
  EqualAreaCylindrical{Datum,latₜₛ}(uconvert(u"m", x), uconvert(u"m", y))
EqualAreaCylindrical{Datum,latₜₛ}(x::Number, y::Number) where {Datum,latₜₛ} =
  EqualAreaCylindrical{Datum,latₜₛ}(addunit(x, u"m"), addunit(y, u"m"))

"""
    Lambert{Datum}(x, y)

Lambert cylindrical equal-area coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
Lambert{WGS84}(1, 1) # add default units
Lambert{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
Lambert{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
Lambert{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54034](https://epsg.io/54034).
"""
const Lambert{Datum} = EqualAreaCylindrical{Datum,0.0u"°"}

typealias(::Type{ESRI{54034}}) = Lambert{WGS84}

"""
    Behrmann{Datum}(x, y)

Behrmann coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
Behrmann{WGS84}(1, 1) # add default units
Behrmann{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
Behrmann{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
Behrmann{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54017](https://epsg.io/54017).
"""
const Behrmann{Datum} = EqualAreaCylindrical{Datum,30.0u"°"}

typealias(::Type{ESRI{54017}}) = Behrmann{WGS84}

"""
    GallPeters{Datum}(x, y)

Gall-Peters coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
GallPeters{WGS84}(1, 1) # add default units
GallPeters{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
GallPeters{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
GallPeters{WGS84}(1.0u"m", 1.0u"m")
```
"""
const GallPeters{Datum} = EqualAreaCylindrical{Datum,45.0u"°"}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EqualAreaCylindrical{Datum,latₜₛ}}, coords::LatLon{Datum}) where {Datum,latₜₛ}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(longitudeₒ(Datum)))
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  l₀ = ustrip(λ₀)
  a = oftype(l, ustrip(majoraxis(🌎)))
  e = oftype(l, eccentricity(🌎))
  e² = oftype(l, eccentricity²(🌎))

  k₀ = cos(ϕₜₛ) / sqrt(1 - e² * sin(ϕₜₛ)^2)
  sinϕ = sin(ϕ)
  esinϕ = e * sinϕ
  q = (1 - e²) * (sinϕ / (1 - esinϕ^2) - (1 / 2e) * log((1 - esinϕ) / (1 + esinϕ)))

  x = a * k₀ * (l - l₀)
  y = a * q / 2k₀

  EqualAreaCylindrical{Datum,latₜₛ}(x * u"m", y * u"m")
end
