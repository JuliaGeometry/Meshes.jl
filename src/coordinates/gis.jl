# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LatLon(lat, lon)

Latitude and longitude in angle units (default to degree).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709)
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

Latitude and longitude in angle units (default to degree)
and altitude in length units (default to meter).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709)
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
