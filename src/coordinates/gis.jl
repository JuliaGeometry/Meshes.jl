# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct EPSG{Code,N,Coords} <: Coordinates{N}
  coords::Coords
end

EPSG{Code,N,Coords}(args...) where {Code,N,Coords} = EPSG{Code,N,Coords}(Coords(args))

Base.isapprox(c₁::T, c₂::T; kwargs...) where {T<:EPSG} =
  all(isapprox(x₁, x₂; kwargs...) for (x₁, x₂) in zip(c₁.coords, c₂.coords))

function Base.show(io::IO, coords::EPSG)
  name = prettyname(coords)
  print(io, "$name(")
  printfields(io, coords.coords, compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", coords::EPSG)
  name = prettyname(coords)
  print(io, "$name coordinates")
  printfields(io, coords.coords)
end

"""
    LatLon(lat, lon)

Latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree).

## Examples

```julia
LatLon(45, 45) # add default units
LatLon(45u"°", 45u"°") # integers are converted converted to floats
LatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
LatLon(45.0u"°", 45.0u"°")
```

See [EPSG:4326](https://epsg.io/4326).
"""
const LatLon{D<:Deg} = EPSG{4326,2,@NamedTuple{lat::D, lon::D}}

LatLon(lat::D, lon::D) where {D<:Deg} = LatLon{float(D)}(lat, lon)
LatLon(lat::Deg, lon::Deg) = LatLon(promote(lat, lon)...)
LatLon(lat::Rad, lon::Rad) = LatLon(rad2deg(lat), rad2deg(lon))
LatLon(lat::Number, lon::Number) = LatLon(addunit(lat, u"°"), addunit(lon, u"°"))

"""
    Mercator(x, y)

Mercator coordinates in length units (default to meter).

## Examples

```julia
Mercator(1, 1) # add default units
Mercator(1u"m", 1u"m") # integers are converted converted to floats
Mercator(1.0u"km", 1.0u"km") # length quantities are converted to meters
Mercator(1.0u"m", 1.0u"m")
```

See [EPSG:3395](https://epsg.io/3395).
"""
const Mercator{M<:Met} = EPSG{3395,2,@NamedTuple{x::M, y::M}}

Mercator(x::M, y::M) where {M<:Met} = Mercator{float(M)}(x, y)
Mercator(x::Met, y::Met) = Mercator(promote(x, y)...)
Mercator(x::Len, y::Len) = Mercator(uconvert(u"m", x), uconvert(u"m", y))
Mercator(x::Number, y::Number) = Mercator(addunit(x, u"m"), addunit(y, u"m"))

"""
    WebMercator(x, y)

WebMercator coordinates in length units (default to meter).

## Examples

```julia
WebMercator(1, 1) # add default units
WebMercator(1u"m", 1u"m") # integers are converted converted to floats
WebMercator(1.0u"km", 1.0u"km") # length quantities are converted to meters
WebMercator(1.0u"m", 1.0u"m")
```

See [EPSG:3857](https://epsg.io/3857).
"""
const WebMercator{M<:Met} = EPSG{3857,2,@NamedTuple{x::M, y::M}}

WebMercator(x::M, y::M) where {M<:Met} = WebMercator{float(M)}(x, y)
WebMercator(x::Met, y::Met) = WebMercator(promote(x, y)...)
WebMercator(x::Len, y::Len) = WebMercator(uconvert(u"m", x), uconvert(u"m", y))
WebMercator(x::Number, y::Number) = WebMercator(addunit(x, u"m"), addunit(y, u"m"))
