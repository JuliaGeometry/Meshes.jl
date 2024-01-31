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
struct LatLon{T<:Deg} <: Coordinates{2}
  lat::T
  lon::T
  LatLon{T}(lat, lon) where {T} = new{float(T)}(lat, lon)
end

LatLon(lat::T, lon::T) where {T<:Deg} = LatLon{T}(lat, lon)
LatLon(lat::Deg, lon::Deg) = LatLon(promote(lat, lon)...)
LatLon(lat::Rad, lon::Rad) = LatLon(rad2deg(lat), rad2deg(lon))
LatLon(lat::Number, lon::Number) = LatLon(lat * u"°", lon * u"°")

"""
    LatLonAlt(lat, lon, alt)

Latitude and longitude in angle units (default to degree)
and altitude in length units (default to meter).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709)
"""
struct LatLonAlt{T<:Deg,A<:Len} <: Coordinates{3}
  lat::T
  lon::T
  alt::A
  LatLonAlt{T,A}(lat, lon, alt) where {T,A} = new{float(T),float(A)}(lat, lon, alt)
end

LatLonAlt(lat::T, lon::T, alt::A) where {T<:Deg,A<:Len} = LatLonAlt{T,A}(lat, lon, alt)
LatLonAlt(lat::Deg, lon::Deg, alt::Len) = LatLonAlt(promote(lat, lon)..., alt)
LatLonAlt(lat::Rad, lon::Rad, alt::Len) = LatLonAlt(rad2deg(lat), rad2deg(lon), alt)
LatLonAlt(lat::Number, lon::Number, alt::Number) = LatLonAlt(lat * u"°", lon * u"°", alt * u"m")

"""
    EastNorth(east, north)

East and north coordinates in length units (default to meter).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
"""
struct EastNorth{T<:Len} <: Coordinates{2}
  east::T
  north::T
  EastNorth{T}(east, north) where {T} = new{float(T)}(east, north)
end

EastNorth(east::T, north::T) where {T<:Len} = EastNorth{T}(east, north)
EastNorth(east::Len, north::Len) = EastNorth(promote(east, north)...)
EastNorth(east::Number, north::Number) = EastNorth(east * u"m", north * u"m")

"""
    WebMercator(x, y)

WebMercator coordinates in length units (default to meter).

## References

* [Web Mercator projection](https://en.wikipedia.org/wiki/Web_Mercator_projection)
"""
struct WebMercator{T<:Quantity} <: Coordinates{2}
  x::T
  y::T
  WebMercator{T}(east, north) where {T} = new{float(T)}(east, north)
end

WebMercator(x::T, y::T) where {T<:Len} = WebMercator{T}(x, y)
WebMercator(x::Len, y::Len) = WebMercator(promote(x, y)...)
WebMercator(x::Number, y::Number) = WebMercator(x * u"m", y * u"m")
