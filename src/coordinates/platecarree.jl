# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
const PlateCarree{M<:Met} = EPSG{32662,2,@NamedTuple{x::M, y::M}}

typealias(::Type{EPSG{32662}}) = PlateCarree

PlateCarree(x::M, y::M) where {M<:Met} = PlateCarree{float(M)}(x, y)
PlateCarree(x::Met, y::Met) = PlateCarree(promote(x, y)...)
PlateCarree(x::Len, y::Len) = PlateCarree(uconvert(u"m", x), uconvert(u"m", y))
PlateCarree(x::Number, y::Number) = PlateCarree(addunit(x, u"m"), addunit(y, u"m"))

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{PlateCarree}, (; coords)::LatLon)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(WGS84.a))
  x = a * l
  y = a * o
  PlateCarree(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, (; coords)::PlateCarree)
  x = coords.x
  y = coords.y
  a = oftype(x, WGS84.a)
  λ = x / a
  ϕ = y / a
  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
