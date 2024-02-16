# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EquidistantCylindrical{Datum,latₜₛ}

Equidistant Cylindrical CRS with a given `Datum` and latitude of true scale `latₜₛ` in degrees.
"""
struct EquidistantCylindrical{Datum,latₜₛ,M<:Met} <: CRS{Datum}
  x::M
  y::M
  EquidistantCylindrical{Datum,latₜₛ}(x::M, y::M) where {Datum,latₜₛ,M<:Met} = new{Datum,latₜₛ,float(M)}(x, y)
end

EquidistantCylindrical{Datum,latₜₛ}(x::Met, y::Met) where {Datum,latₜₛ} =
  EquidistantCylindrical{Datum,latₜₛ}(promote(x, y)...)
EquidistantCylindrical{Datum,latₜₛ}(x::Len, y::Len) where {Datum,latₜₛ} =
  EquidistantCylindrical{Datum,latₜₛ}(uconvert(u"m", x), uconvert(u"m", y))
EquidistantCylindrical{Datum,latₜₛ}(x::Number, y::Number) where {Datum,latₜₛ} =
  EquidistantCylindrical{Datum,latₜₛ}(addunit(x, u"m"), addunit(y, u"m"))

"""
    PlateCarree{Datum}(x, y)

Plate Carrée coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
PlateCarree{WGS84}(1, 1) # add default units
PlateCarree{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
PlateCarree{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
PlateCarree{WGS84}(1.0u"m", 1.0u"m")
```

See [EPSG:32662](https://epsg.io/32662).
"""
const PlateCarree{Datum} = EquidistantCylindrical{Datum,0.0u"°"}

typealias(::Type{EPSG{32662}}) = PlateCarree{WGS84}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EquidistantCylindrical{Datum,latₜₛ}}, coords::LatLon{Datum}) where {Datum,latₜₛ}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  x = a * l * cos(ϕₜₛ)
  y = a * o

  EquidistantCylindrical{Datum,latₜₛ}(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon{Datum}}, coords::EquidistantCylindrical{Datum,latₜₛ}) where {Datum,latₜₛ}
  🌎 = ellipsoid(Datum)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(🌎))
  ϕₜₛ = numconvert(numtype(x), deg2rad(latₜₛ))

  λ = x / (cos(ϕₜₛ) * a)
  ϕ = y / a

  LatLon{Datum}(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
