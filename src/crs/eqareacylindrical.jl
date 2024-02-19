# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EqualAreaCylindrical{latₜₛ,Datum}

Equal Area Cylindrical CRS with latitude of true scale `latₜₛ` in degrees and a given `Datum`.
"""
struct EqualAreaCylindrical{latₜₛ,Datum,M<:Met} <: CRS{Datum}
  x::M
  y::M
  EqualAreaCylindrical{latₜₛ,Datum}(x::M, y::M) where {latₜₛ,Datum,M<:Met} = new{latₜₛ,Datum,float(M)}(x, y)
end

EqualAreaCylindrical{latₜₛ,Datum}(x::Met, y::Met) where {latₜₛ,Datum} =
  EqualAreaCylindrical{latₜₛ,Datum}(promote(x, y)...)
EqualAreaCylindrical{latₜₛ,Datum}(x::Len, y::Len) where {latₜₛ,Datum} =
  EqualAreaCylindrical{latₜₛ,Datum}(uconvert(u"m", x), uconvert(u"m", y))
EqualAreaCylindrical{latₜₛ,Datum}(x::Number, y::Number) where {latₜₛ,Datum} =
  EqualAreaCylindrical{latₜₛ,Datum}(addunit(x, u"m"), addunit(y, u"m"))

EqualAreaCylindrical{latₜₛ}(args...) where {latₜₛ} = EqualAreaCylindrical{latₜₛ,WGS84}(args...)

"""
    Lambert(x, y)
    Lambert{Datum}(x, y)

Lambert cylindrical equal-area coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
Lambert(1, 1) # add default units
Lambert(1u"m", 1u"m") # integers are converted converted to floats
Lambert(1.0u"km", 1.0u"km") # length quantities are converted to meters
Lambert(1.0u"m", 1.0u"m")
Lambert{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54034](https://epsg.io/54034).
"""
const Lambert{Datum} = EqualAreaCylindrical{0.0u"°",Datum}

"""
    Behrmann(x, y)
    Behrmann{Datum}(x, y)

Behrmann coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
Behrmann(1, 1) # add default units
Behrmann(1u"m", 1u"m") # integers are converted converted to floats
Behrmann(1.0u"km", 1.0u"km") # length quantities are converted to meters
Behrmann(1.0u"m", 1.0u"m")
Behrmann{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54017](https://epsg.io/54017).
"""
const Behrmann{Datum} = EqualAreaCylindrical{30.0u"°",Datum}

"""
    GallPeters(x, y)
    GallPeters{Datum}(x, y)

Gall-Peters coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
GallPeters(1, 1) # add default units
GallPeters(1u"m", 1u"m") # integers are converted converted to floats
GallPeters(1.0u"km", 1.0u"km") # length quantities are converted to meters
GallPeters(1.0u"m", 1.0u"m")
GallPeters{WGS84}(1.0u"m", 1.0u"m")
```
"""
const GallPeters{Datum} = EqualAreaCylindrical{45.0u"°",Datum}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EqualAreaCylindrical{latₜₛ,Datum}}, coords::LatLon{Geodetic,Datum}) where {latₜₛ,Datum}
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

  EqualAreaCylindrical{latₜₛ,Datum}(x * u"m", y * u"m")
end
