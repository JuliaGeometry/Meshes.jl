# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LatLon(lat, lon)

Latitude and longitude in angular units (default to degree).

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
struct LatLon{D<:Deg} <: Coordinates{2}
  lat::D
  lon::D
  LatLon{D}(lat, lon) where {D} = new{float(D)}(lat, lon)
end

LatLon(lat::D, lon::D) where {D<:Deg} = LatLon{D}(lat, lon)
LatLon(lat::Deg, lon::Deg) = LatLon(promote(lat, lon)...)
LatLon(lat::Rad, lon::Rad) = LatLon(rad2deg(lat), rad2deg(lon))
LatLon(lat::Number, lon::Number) = LatLon(addunit(lat, u"°"), addunit(lon, u"°"))

"""
    LatLonAlt(lat, lon, alt)

Latitude and longitude in angular units (default to degree)
and altitude in length units (default to meter).

## Examples

```julia
LatLonAlt(45, 45, 1) # add default units
LatLonAlt(45u"°", 45u"°", 1u"m") # integers are converted converted to floats
LatLonAlt((π/4)u"rad", (π/4)u"rad", 1.0u"m") # radians are converted to degrees
LatLonAlt(45.0u"°", 45.0u"°", 1.0u"km")
```

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709:2022](https://www.iso.org/standard/75147.html)
"""
struct LatLonAlt{D<:Deg,L<:Len} <: Coordinates{3}
  lat::D
  lon::D
  alt::L
  LatLonAlt{D,L}(lat, lon, alt) where {D,L} = new{float(D),float(L)}(lat, lon, alt)
end

LatLonAlt(lat::D, lon::D, alt::L) where {D<:Deg,L<:Len} = LatLonAlt{D,L}(lat, lon, alt)
LatLonAlt(lat::Deg, lon::Deg, alt::Len) = LatLonAlt(promote(lat, lon)..., alt)
LatLonAlt(lat::Rad, lon::Rad, alt::Len) = LatLonAlt(rad2deg(lat), rad2deg(lon), alt)
LatLonAlt(lat::Number, lon::Number, alt::Number) = LatLonAlt(addunit(lat, u"°"), addunit(lon, u"°"), addunit(alt, u"m"))

"""
    EastNorth(east, north)

East and north coordinates in length units (default to meter).

## Examples

```julia
EastNorth(1, 1) # add default units
EastNorth(1"m", 1"m") # integers are converted converted to floats
EastNorth(1.0"km", 1.0"km")
```

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
"""
struct EastNorth{L<:Len} <: Coordinates{2}
  east::L
  north::L
  EastNorth{L}(east, north) where {L} = new{float(L)}(east, north)
end

EastNorth(east::L, north::L) where {L<:Len} = EastNorth{L}(east, north)
EastNorth(east::Len, north::Len) = EastNorth(promote(east, north)...)
EastNorth(east::Number, north::Number) = EastNorth(addunit(east, u"m"), addunit(north, u"m"))

"""
    WebMercator(x, y)

WebMercator coordinates in length units (default to meter).

## Examples

```julia
WebMercator(1, 1) # add default units
WebMercator(1"m", 1"m") # integers are converted converted to floats
WebMercator(1.0"km", 1.0"km")
```

## References

* [Web Mercator projection](https://en.wikipedia.org/wiki/Web_Mercator_projection)
"""
struct WebMercator{L<:Len} <: Coordinates{2}
  x::L
  y::L
  WebMercator{L}(east, north) where {L} = new{float(L)}(east, north)
end

WebMercator(x::L, y::L) where {L<:Len} = WebMercator{L}(x, y)
WebMercator(x::Len, y::Len) = WebMercator(promote(x, y)...)
WebMercator(x::Number, y::Number) = WebMercator(addunit(x, u"m"), addunit(y, u"m"))
