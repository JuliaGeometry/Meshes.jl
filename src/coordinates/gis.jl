# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LatLon(lat, lon)

Latitude and longitude in degrees.

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709)
"""
struct LatLon{T<:Quantity} <: Coordinates{2,T}
  lat::T
  lon::T
  function LatLon{T}(lat, lon) where {T<:Quantity}
    if unit(T) ≠ u"°"
      throw(ArgumentError("the units of `lat` and `lon` must be degrees"))
    end
    new{float(T)}(lat, lon)
  end
end

LatLon(lat::T, lon::T) where {T<:Quantity} = LatLon{T}(lat, lon)
LatLon(lat::Quantity, lon::Quantity) = LatLon(promote(lat, lon)...)
LatLon(lat::Number, lon::Number) = LatLon(lat * u"°", lon * u"°")

"""
    LatLonAlt(lat, lon, alt)

Latitude and longitude in degrees and altitude in length units (default to meter).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
* [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709)
"""
struct LatLonAlt{T<:Quantity,A<:Quantity} <: Coordinates{3,T}
  lat::T
  lon::T
  alt::A
  function LatLonAlt{T,A}(lat, lon, alt) where {T<:Quantity,A<:Quantity}
    if unit(T) ≠ u"°"
      throw(ArgumentError("the units of `lat` and `lon` must be degrees"))
    end
    if dimension(A) ≠ u"𝐋"
      throw(ArgumentError("the unit of `alt` must be a length unit"))
    end
    new{float(T),float(A)}(lat, lon, alt)
  end
end

LatLonAlt(lat::T, lon::T, alt::A) where {T<:Quantity,A<:Quantity} = LatLonAlt{T,A}(lat, lon, alt)
LatLonAlt(lat::Quantity, lon::Quantity, alt::Quantity) = LatLonAlt(promote(lat, lon)..., alt)
LatLonAlt(lat::Number, lon::Number, alt::Number) = LatLonAlt(lat * u"°", lon * u"°", alt * u"m")

"""
    EastNorth(east, north)

East and north coordinates in length units (default to meter).

## References

* [Geographic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system)
"""
struct EastNorth{T<:Quantity} <: Coordinates{2,T}
  east::T
  north::T
  function EastNorth{T}(east, north) where {T<:Quantity}
    if dimension(T) ≠ u"𝐋"
      throw(ArgumentError("the units of `east` and `north` must be length units"))
    end
    new{float(T)}(east, north)
  end
end

EastNorth(east::T, north::T) where {T<:Quantity} = EastNorth{T}(east, north)
EastNorth(east::Quantity, north::Quantity) = EastNorth(promote(east, north)...)
EastNorth(east::Number, north::Number) = EastNorth(east * u"m", north * u"m")

"""
    WebMercator(x, y)

WebMercator coordinates in length units (default to meter).

## References

* [Web Mercator projection](https://en.wikipedia.org/wiki/Web_Mercator_projection)
"""
struct WebMercator{T<:Quantity} <: Coordinates{2,T}
  x::T
  y::T
  function WebMercator{T}(x, y) where {T<:Quantity}
    if dimension(T) ≠ u"𝐋"
      throw(ArgumentError("the units of `x` and `y` must be length units"))
    end
    new{float(T)}(x, y)
  end
end

WebMercator(x::T, y::T) where {T<:Quantity} = WebMercator{T}(x, y)
WebMercator(x::Quantity, y::Quantity) = WebMercator(promote(x, y)...)
WebMercator(x::Number, y::Number) = WebMercator(x * u"m", y * u"m")
