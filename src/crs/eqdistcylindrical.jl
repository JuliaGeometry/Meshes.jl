# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EquidistantCylindrical{latₜₛ}

Equidistant Cylindrical CRS with latitude of true scale `latₜₛ` in degrees.
"""
const EquidistantCylindrical{latₜₛ,M<:Met} = CRS{:EDC,@NamedTuple{x::M, y::M},WGS84,latₜₛ}

EquidistantCylindrical{latₜₛ}(x::M, y::M) where {latₜₛ,M<:Met} = EquidistantCylindrical{latₜₛ,float(M)}(x, y)
EquidistantCylindrical{latₜₛ}(x::Met, y::Met) where {latₜₛ} = EquidistantCylindrical{latₜₛ}(promote(x, y)...)
EquidistantCylindrical{latₜₛ}(x::Len, y::Len) where {latₜₛ} =
  EquidistantCylindrical{latₜₛ}(uconvert(u"m", x), uconvert(u"m", y))
EquidistantCylindrical{latₜₛ}(x::Number, y::Number) where {latₜₛ} =
  EquidistantCylindrical{latₜₛ}(addunit(x, u"m"), addunit(y, u"m"))

"""
    PlateCarree(x, y)

Plate Carrée coordinates in length units (default to meter).

## Examples

```julia
PlateCarree(1, 1) # add default units
PlateCarree(1u"m", 1u"m") # integers are converted converted to floats
PlateCarree(1.0u"km", 1.0u"km") # length quantities are converted to meters
PlateCarree(1.0u"m", 1.0u"m")
```

See [EPSG:32662](https://epsg.io/32662).
"""
const PlateCarree = EquidistantCylindrical{0.0u"°"}

typealias(::Type{EPSG{32662}}) = PlateCarree

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EquidistantCylindrical{latₜₛ}}, coords::LatLon) where {latₜₛ}
  🌎 = ellipsoid(EquidistantCylindrical{latₜₛ})
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  x = a * l * cos(ϕₜₛ)
  y = a * o

  EquidistantCylindrical{latₜₛ}(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, coords::EquidistantCylindrical{latₜₛ}) where {latₜₛ}
  🌎 = ellipsoid(coords)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(🌎))
  ϕₜₛ = numconvert(numtype(x), deg2rad(latₜₛ))

  λ = x / (cos(ϕₜₛ) * a)
  ϕ = y / a

  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
