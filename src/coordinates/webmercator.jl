# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WebMercator(x, y)

Web Mercator coordinates in length units (default to meter).

## Examples

```julia
WebMercator(1, 1) # add default units
WebMercator(1u"m", 1u"m") # integers are converted converted to floats
WebMercator(1.0u"km", 1.0u"km") # length quantities are converted to meters
WebMercator(1.0u"m", 1.0u"m")
```

See [EPSG:3857](https://epsg.io/3857).
"""
const WebMercator{M<:Met} = CRS{EPSG{3857},@NamedTuple{x::M, y::M},WGS84}

typealias(::Type{EPSG{3857}}) = WebMercator

WebMercator(x::M, y::M) where {M<:Met} = WebMercator{float(M)}(x, y)
WebMercator(x::Met, y::Met) = WebMercator(promote(x, y)...)
WebMercator(x::Len, y::Len) = WebMercator(uconvert(u"m", x), uconvert(u"m", y))
WebMercator(x::Number, y::Number) = WebMercator(addunit(x, u"m"), addunit(y, u"m"))

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{WebMercator}, coords::LatLon)
  ellip = ellipsoid(coords)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(ellip)))
  x = a * l
  y = a * asinh(tan(ϕ))
  WebMercator(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, coords::WebMercator)
  ellip = ellipsoid(coords)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(ellip))
  λ = x / a
  ϕ = atan(sinh(y / a))
  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
