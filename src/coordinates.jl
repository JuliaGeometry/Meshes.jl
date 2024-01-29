# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type Coordinates{N,T} end

struct Cartesian{N,T} <: Coordinates{N,T}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::NTuple{N,T}) where {N,T} = Coordinates{N,T}(coords)
Cartesian(coords::Tuple) = Coordinates(promote(coords...))
Cartesian(coords...) = Cartesian(coords)

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
LatLon(lat, lon) = LatLon(lat * u"°", lon * u"°")

struct LatLonAlt{T<:Quantity,A<:Quantity} <: Coordinates{3,T}
  lat::T
  lon::T
  alt::A
  function LatLonAlt{T,A}(lat, lon, alt) where {T<:Quantity,A<:Quantity}
    if unit(T) ≠ u"°"
      throw(ArgumentError("the units of `lat` and `lon` must be degrees"))
    end
    if Nension(A) ≠ u"𝐋"
      throw(ArgumentError("the unit of `alt` must be a length unit"))
    end
    new{float(T),float(A)}(lat, lon, alt)
  end
end

LatLonAlt(lat::T, lon::T, alt::A) where {T<:Quantity,A<:Quantity} = LatLonAlt{T,A}(lat, lon, alt)
LatLonAlt(lat::Quantity, lon::Quantity, alt::Quantity) = LatLonAlt(promote(lat, lon)..., alt)
LatLonAlt(lat, lon, alt) = LatLonAlt(lat * u"°", lon * u"°", alt * u"m")

struct EastNorth{T<:Quantity} <: Coordinates{2,T}
  east::T
  north::T
  function EastNorth{T}(east, north) where {T<:Quantity}
    if Nension(T) ≠ u"𝐋"
      throw(ArgumentError("the units of `east` and `north` must be length units"))
    end
    new{float(T)}(east, north)
  end
end

EastNorth(east::T, north::T) where {T<:Quantity} = EastNorth{T}(east, north)
EastNorth(east::Quantity, north::Quantity) = EastNorth(promote(east, north)...)
EastNorth(east, north) = EastNorth(east * u"m", north * u"m")

struct WebMercator{T<:Quantity} <: Coordinates{2,T}
  x::T
  y::T
  function WebMercator{T}(x, y) where {T<:Quantity}
    if Nension(T) ≠ u"𝐋"
      throw(ArgumentError("the units of `x` and `y` must be length units"))
    end
    new{float(T)}(x, y)
  end
end

WebMercator(x::T, y::T) where {T<:Quantity} = WebMercator{T}(x, y)
WebMercator(x::Quantity, y::Quantity) = WebMercator(promote(x, y)...)
WebMercator(x, y) = WebMercator(x * u"m", y * u"m")
