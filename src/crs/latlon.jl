# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type LatitudeLongitude{Datum} <: CRS{Datum} end

"""
    LatLon(lat, lon)
    LatLon{Datum}(lat, lon)
    GeodeticLatLon(lat, lon)
    GeodeticLatLon{Datum}(lat, lon)

Geodetic latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`).

`LatLon` is an alias to `GeodeticLatLon`.

## Examples

```julia
LatLon(45, 45) # add default units
LatLon(45u"°", 45u"°") # integers are converted converted to floats
LatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
LatLon(45.0u"°", 45.0u"°")
LatLon{WGS84}(45.0u"°", 45.0u"°")
```

See [EPSG:4326](https://epsg.io/4326).
"""
struct GeodeticLatLon{Datum,D<:Deg} <: LatitudeLongitude{Datum}
  lat::D
  lon::D
  GeodeticLatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

GeodeticLatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = GeodeticLatLon{Datum}(promote(lat, lon)...)
GeodeticLatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = GeodeticLatLon{Datum}(rad2deg(lat), rad2deg(lon))
GeodeticLatLon{Datum}(lat::Number, lon::Number) where {Datum} =
  GeodeticLatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

GeodeticLatLon(args...) = GeodeticLatLon{WGS84}(args...)

const LatLon = GeodeticLatLon

"""
    GeocentricLatLon(lat, lon)
    GeocentricLatLon{Datum}(lat, lon)

Geocentric latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
GeocentricLatLon(45, 45) # add default units
GeocentricLatLon(45u"°", 45u"°") # integers are converted converted to floats
GeocentricLatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
GeocentricLatLon(45.0u"°", 45.0u"°")
GeocentricLatLon{WGS84}(45.0u"°", 45.0u"°")
```
"""
struct GeocentricLatLon{Datum,D<:Deg} <: LatitudeLongitude{Datum}
  lat::D
  lon::D
  GeocentricLatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

GeocentricLatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = GeocentricLatLon{Datum}(promote(lat, lon)...)
GeocentricLatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = GeocentricLatLon{Datum}(rad2deg(lat), rad2deg(lon))
GeocentricLatLon{Datum}(lat::Number, lon::Number) where {Datum} =
  GeocentricLatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

GeocentricLatLon(args...) = GeocentricLatLon{WGS84}(args...)

"""
    AuthalicLatLon(lat, lon)
    AuthalicLatLon{Datum}(lat, lon)

Authalic latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
AuthalicLatLon(45, 45) # add default units
AuthalicLatLon(45u"°", 45u"°") # integers are converted converted to floats
AuthalicLatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
AuthalicLatLon(45.0u"°", 45.0u"°")
AuthalicLatLon{WGS84}(45.0u"°", 45.0u"°")
```
"""
struct AuthalicLatLon{Datum,D<:Deg} <: LatitudeLongitude{Datum}
  lat::D
  lon::D
  AuthalicLatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

AuthalicLatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = AuthalicLatLon{Datum}(promote(lat, lon)...)
AuthalicLatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = AuthalicLatLon{Datum}(rad2deg(lat), rad2deg(lon))
AuthalicLatLon{Datum}(lat::Number, lon::Number) where {Datum} =
  AuthalicLatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

AuthalicLatLon(args...) = AuthalicLatLon{WGS84}(args...)

# ------------
# CONVERSIONS
# ------------

# Adapted from PROJ coordinate transformation software
# Initial PROJ 4.3 public domain code was put as Frank Warmerdam as copyright
# holder, but he didn't mean to imply he did the work. Essentially all work was
# done by Gerald Evenden.

# reference code: https://github.com/OSGeo/PROJ/blob/master/src/4D_api.cpp#L774

function Base.convert(::Type{GeocentricLatLon{Datum}}, (; lat, lon)::LatLon{Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan((1 - e²) * tan(lat)))
  GeocentricLatLon{Datum}(lat′ * u"°", lon)
end

function Base.convert(::Type{LatLon{Datum}}, (; lat, lon)::GeocentricLatLon{Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan(1 / (1 - e²) * tan(lat)))
  LatLon{Datum}(lat′ * u"°", lon)
end

# reference code: https://github.com/OSGeo/PROJ/blob/master/src/projections/healpix.cpp#L230
# reference formula: https://mathworld.wolfram.com/AuthalicLatitude.html

function Base.convert(::Type{AuthalicLatLon{Datum}}, (; lat, lon)::LatLon{Datum}) where {Datum}
  🌎 = ellipsoid(Datum)
  l = ustrip(lat)
  e = oftype(l, eccentricity(🌎))
  e² = oftype(l, eccentricity²(🌎))
  ome² = 1 - e²
  sinϕ = sin(lat)
  esinϕ = e * sinϕ

  q = ome² * (sinϕ / (1 - esinϕ^2) - (1 / 2e) * log((1 - esinϕ) / (1 + esinϕ)))
  # same formula as q, but ϕ = 90°
  qₚ = ome² * (1 / ome² - (1 / 2e) * log((1 - e) / (1 + e)))
  qqₚ⁻¹ = q / qₚ

  if abs(qqₚ⁻¹) > 1
    qqₚ⁻¹ = sign(qqₚ⁻¹) # rounding error
  end

  lat′ = rad2deg(asin(qqₚ⁻¹))
  AuthalicLatLon{Datum}(lat′ * u"°", lon)
end

# reference code: https://github.com/OSGeo/PROJ/blob/master/src/auth.cpp
# reference formula: https://mathworld.wolfram.com/AuthalicLatitude.html

const _P₁₁ = 0.33333333333333333333 # 1 / 3
const _P₁₂ = 0.17222222222222222222 # 31 / 180
const _P₁₃ = 0.10257936507936507937 # 517 / 5040
const _P₂₁ = 0.06388888888888888888 # 23 / 360
const _P₂₂ = 0.06640211640211640212 # 251 / 3780
const _P₃₁ = 0.01677689594356261023 # 761 / 45360

function authlat(β, e²)
  e⁴ = e²^2
  e⁶ = e²^3
  P₁₁ = oftype(β, _P₁₁)
  P₁₂ = oftype(β, _P₁₂)
  P₁₃ = oftype(β, _P₁₃)
  P₂₁ = oftype(β, _P₂₁)
  P₂₂ = oftype(β, _P₂₂)
  P₃₁ = oftype(β, _P₃₁)
  β + (P₁₁ * e² + P₁₂ * e⁴ + P₁₃ * e⁶) * sin(2β) + (P₂₁ * e⁴ + P₂₂ * e⁶) * sin(4β) + (P₃₁ * e⁶) * sin(6β)
end

function Base.convert(::Type{LatLon{Datum}}, (; lat, lon)::AuthalicLatLon{Datum}) where {Datum}
  l = ustrip(deg2rad(lat))
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(authlat(l, e²))
  LatLon{Datum}(lat′ * u"°", lon)
end
