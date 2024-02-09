# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EquidistantCylindrical{ID,latₜₛ}

Equidistant Cylindrical CRS with identifier `ID` and latitude of true scale `latₜₛ` in degrees.
"""
const EquidistantCylindrical{ID,latₜₛ,M<:Met} = CRS{ID,@NamedTuple{x::M, y::M},WGS84,latₜₛ}

EquidistantCylindrical{ID,latₜₛ}(x::M, y::M) where {ID,latₜₛ,M<:Met} = EquidistantCylindrical{ID,latₜₛ,float(M)}(x, y)
EquidistantCylindrical{ID,latₜₛ}(x::Met, y::Met) where {ID,latₜₛ} = EquidistantCylindrical{ID,latₜₛ}(promote(x, y)...)
EquidistantCylindrical{ID,latₜₛ}(x::Len, y::Len) where {ID,latₜₛ} =
  EquidistantCylindrical{ID,latₜₛ}(uconvert(u"m", x), uconvert(u"m", y))
EquidistantCylindrical{ID,latₜₛ}(x::Number, y::Number) where {ID,latₜₛ} =
  EquidistantCylindrical{ID,latₜₛ}(addunit(x, u"m"), addunit(y, u"m"))

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
const PlateCarree = EquidistantCylindrical{EPSG{32662},0.0u"°"}

typealias(::Type{EPSG{32662}}) = PlateCarree

function Base.convert(::Type{EquidistantCylindrical{ID,latₜₛ}}, coords::LatLon) where {ID,latₜₛ}
  🌎 = ellipsoid(coords)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  x = a * l * cos(ϕₜₛ)
  y = a * o

  EquidistantCylindrical{ID,latₜₛ}(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, coords::EquidistantCylindrical{ID,latₜₛ}) where {ID,latₜₛ}
  🌎 = ellipsoid(coords)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(🌎))
  ϕₜₛ = numconvert(numtype(x), deg2rad(latₜₛ))

  λ = x / (cos(ϕₜₛ) * a)
  ϕ = y / a

  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
