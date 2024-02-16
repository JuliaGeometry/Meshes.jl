# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Orthographic{Datum,lat₀,lon₀,S}

Orthographic CRS with with a given `Datum`, latitude origin `lat₀` and longitude origin `lon₀` in degrees
and spherical mode `S` enabled or not.
"""
struct Orthographic{Datum,lat₀,lon₀,S,M<:Met} <: CRS{Datum}
  x::M
  y::M
  Orthographic{Datum,lat₀,lon₀,S}(x::M, y::M) where {Datum,lat₀,lon₀,S,M<:Met} = new{Datum,lat₀,lon₀,S,float(M)}(x, y)
end

Orthographic{Datum,lat₀,lon₀,S}(x::Met, y::Met) where {Datum,lat₀,lon₀,S} =
  Orthographic{Datum,lat₀,lon₀,S}(promote(x, y)...)
Orthographic{Datum,lat₀,lon₀,S}(x::Len, y::Len) where {Datum,lat₀,lon₀,S} =
  Orthographic{Datum,lat₀,lon₀,S}(uconvert(u"m", x), uconvert(u"m", y))
Orthographic{Datum,lat₀,lon₀,S}(x::Number, y::Number) where {Datum,lat₀,lon₀,S} =
  Orthographic{Datum,lat₀,lon₀,S}(addunit(x, u"m"), addunit(y, u"m"))

"""
    OrthoNorth{Datum}(x, y)

Orthographic North Pole coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
OrthoNorth{WGS84}(1, 1) # add default units
OrthoNorth{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
OrthoNorth{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
OrthoNorth{WGS84}(1.0u"m", 1.0u"m")
```
"""
const OrthoNorth{Datum} = Orthographic{Datum,90.0u"°",0.0u"°",false}

"""
    OrthoSouth{Datum}(x, y)

Orthographic South Pole coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
OrthoSouth{WGS84}(1, 1) # add default units
OrthoSouth{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
OrthoSouth{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
OrthoSouth{WGS84}(1.0u"m", 1.0u"m")
```
"""
const OrthoSouth{Datum} = Orthographic{Datum,-90.0u"°",0.0u"°",false}

typealias(::Type{ESRI{102035}}) = Orthographic{WGS84,90.0u"°",0.0u"°",true}

typealias(::Type{ESRI{102037}}) = Orthographic{WGS84,-90.0u"°",0.0u"°",true}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Orthographic{Datum,lat₀,lon₀,false}}, coords::LatLon{Datum}) where {Datum,lat₀,lon₀}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(lon₀))
  ϕ₀ = oftype(ϕ, deg2rad(lat₀))
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  e² = oftype(l, eccentricity²(🌎))

  sinϕ = sin(ϕ)
  cosϕ = cos(ϕ)
  sinϕ₀ = sin(ϕ₀)
  cosϕ₀ = cos(ϕ₀)
  ν = 1 / sqrt(1 - e² * sinϕ^2)
  ν₀ = 1 / sqrt(1 - e² * sinϕ₀^2)

  x = a * (ν * cosϕ * sin(λ - λ₀))
  y = a * (ν * (sinϕ * cosϕ₀ - cosϕ * sinϕ₀ * cos(λ - λ₀)) + e² * (ν₀ * sinϕ₀ - ν * sinϕ) * cosϕ₀)

  Orthographic{Datum,lat₀,lon₀,false}(x * u"m", y * u"m")
end

function Base.convert(::Type{Orthographic{Datum,lat₀,lon₀,true}}, coords::LatLon{Datum}) where {Datum,lat₀,lon₀}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(lon₀))
  ϕ₀ = oftype(ϕ, deg2rad(lat₀))
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  cosϕ = cos(ϕ)
  x = a * cosϕ * sin(λ - λ₀)
  y = a * (sin(ϕ) * cos(ϕ₀) - cosϕ * sin(ϕ₀) * cos(λ - λ₀))

  Orthographic{Datum,lat₀,lon₀,true}(x * u"m", y * u"m")
end
