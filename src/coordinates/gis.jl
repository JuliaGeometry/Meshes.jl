# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct EPSG{N,Code,Coords} <: Coordinates{N}
  coords::Coords
end

EPSG{N,Code,Coords}(; kwargs...) where {N,Code,Coords} = EPSG{N,Code,Coords}(values(kwargs))

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

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709:2022](https://www.iso.org/standard/75147.html)
"""
const LatLon{D<:Deg} = EPSG{2,4326,@NamedTuple{lat::D, lon::D}}

LatLon(lat::D, lon::D) where {D<:Deg} = LatLon{float(D)}(; lat, lon)
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
Mercator(1.0u"km", 1.0u"km")
```

## References

* [Mercator projection](https://en.wikipedia.org/wiki/Mercator_projection)
"""
const Mercator{L<:Len} = EPSG{2,3395,@NamedTuple{x::L, y::L}}

Mercator(x::L, y::L) where {L<:Len} = Mercator{float(L)}(; x, y)
Mercator(x::Len, y::Len) = Mercator(promote(x, y)...)
Mercator(x::Number, y::Number) = Mercator(addunit(x, u"m"), addunit(y, u"m"))

"""
    WebMercator(x, y)

WebMercator coordinates in length units (default to meter).

## Examples

```julia
WebMercator(1, 1) # add default units
WebMercator(1u"m", 1u"m") # integers are converted converted to floats
WebMercator(1.0u"km", 1.0u"km")
```

## References

* [Web Mercator projection](https://en.wikipedia.org/wiki/Web_Mercator_projection)
"""
const WebMercator{L<:Len} = EPSG{2,3857,@NamedTuple{x::L, y::L}}

WebMercator(x::L, y::L) where {L<:Len} = WebMercator{float(L)}(; x, y)
WebMercator(x::Len, y::Len) = WebMercator(promote(x, y)...)
WebMercator(x::Number, y::Number) = WebMercator(addunit(x, u"m"), addunit(y, u"m"))
