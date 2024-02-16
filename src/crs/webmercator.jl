# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WebMercator{Datum}(x, y)

Web Mercator coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
WebMercator{WGS84}(1, 1) # add default units
WebMercator{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
WebMercator{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
WebMercator{WGS84}(1.0u"m", 1.0u"m")
```

See [EPSG:3857](https://epsg.io/3857).
"""
struct WebMercator{Datum,M<:Met} <: CRS{Datum}
  x::M
  y::M
  WebMercator{Datum}(x::M, y::M) where {Datum,M<:Met} = new{Datum,float(M)}(x, y)
end

typealias(::Type{EPSG{3857}}) = WebMercator{WGS84}

WebMercator{Datum}(x::Met, y::Met) where {Datum} = WebMercator{Datum}(promote(x, y)...)
WebMercator{Datum}(x::Len, y::Len) where {Datum} = WebMercator{Datum}(uconvert(u"m", x), uconvert(u"m", y))
WebMercator{Datum}(x::Number, y::Number) where {Datum} = WebMercator{Datum}(addunit(x, u"m"), addunit(y, u"m"))

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{WebMercator{Datum}}, coords::LatLon{Datum}) where {Datum}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  x = a * l
  y = a * asinh(tan(ϕ))
  WebMercator{Datum}(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon{Datum}}, coords::WebMercator{Datum}) where {Datum}
  🌎 = ellipsoid(Datum)
  x = coords.x
  y = coords.y
  a = oftype(x, majoraxis(🌎))
  λ = x / a
  ϕ = atan(sinh(y / a))
  LatLon{Datum}(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end
