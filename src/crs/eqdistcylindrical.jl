# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EquidistantCylindrical{latₜₛ,Datum}

Equidistant Cylindrical CRS with latitude of true scale `latₜₛ` in degrees and a given `Datum`.
"""
struct EquidistantCylindrical{latₜₛ,Datum,M<:Met} <: CRS{Datum}
  x::M
  y::M
  EquidistantCylindrical{latₜₛ,Datum}(x::M, y::M) where {latₜₛ,Datum,M<:Met} = new{latₜₛ,Datum,float(M)}(x, y)
end

EquidistantCylindrical{latₜₛ,Datum}(x::Met, y::Met) where {latₜₛ,Datum} =
  EquidistantCylindrical{latₜₛ,Datum}(promote(x, y)...)
EquidistantCylindrical{latₜₛ,Datum}(x::Len, y::Len) where {latₜₛ,Datum} =
  EquidistantCylindrical{latₜₛ,Datum}(uconvert(u"m", x), uconvert(u"m", y))
EquidistantCylindrical{latₜₛ,Datum}(x::Number, y::Number) where {latₜₛ,Datum} =
  EquidistantCylindrical{latₜₛ,Datum}(addunit(x, u"m"), addunit(y, u"m"))

EquidistantCylindrical{latₜₛ}(args...) where {latₜₛ} = EquidistantCylindrical{latₜₛ,WGS84}(args...)

"""
    PlateCarree(x, y)
    PlateCarree{Datum}(x, y)

Plate Carrée coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
PlateCarree(1, 1) # add default units
PlateCarree(1u"m", 1u"m") # integers are converted converted to floats
PlateCarree(1.0u"km", 1.0u"km") # length quantities are converted to meters
PlateCarree(1.0u"m", 1.0u"m")
PlateCarree{WGS84}(1.0u"m", 1.0u"m")
```

See [EPSG:32662](https://epsg.io/32662).
"""
const PlateCarree{Datum} = EquidistantCylindrical{0.0u"°",Datum}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EquidistantCylindrical{latₜₛ,Datum}}, coords::LatLon{Datum}) where {latₜₛ,Datum}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  x = a * l * cos(ϕₜₛ)
  y = a * o

  EquidistantCylindrical{latₜₛ,Datum}(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon{Datum}}, coords::EquidistantCylindrical{latₜₛ,Datum}) where {latₜₛ,Datum}
  🌎 = ellipsoid(Datum)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(🌎))
  ϕₜₛ = numconvert(numtype(x), deg2rad(latₜₛ))

  λ = x / (cos(ϕₜₛ) * a)
  ϕ = y / a

  LatLon{Datum}(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
