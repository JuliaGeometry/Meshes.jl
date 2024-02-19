abstract type LatitudeLongitude{Datum} <: CRS{Datum} end

"""
    LatLon(lat, lon)
    LatLon{Datum}(lat, lon)
    GeodeticLatLon(lat, lon)
    GeodeticLatLon{Datum}(lat, lon)

Geodetic latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`). `LatLon` is an alias to `GeocentricLatLon` and the recommended constructor.

## Examples

```julia
LatLon(45, 45) # add default units
LatLon(45u"°", 45u"°") # integers are converted converted to floats
LatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
LatLon(45.0u"°", 45.0u"°")
LatLon{WGS84}(45.0u"°", 45.0u"°")
GeodeticLatLon(45.0u"°", 45.0u"°")
GeodeticLatLon{WGS84}(45.0u"°", 45.0u"°")
```

See [EPSG:4326](https://epsg.io/4326).
"""
struct GeodeticLatLon{Datum,D<:Deg} <: LatitudeLongitude{Datum}
  lat::D
  lon::D
  GeodeticLatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

GeodeticLatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = GeodeticLatLon{Datum}(promote(lat, lon)...)
GeodeticLatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = GeodeticLatLon{Datum}(rad2deg(lat), rad2deg(lon))
GeodeticLatLon{Datum}(lat::Number, lon::Number) where {Datum} =
  GeodeticLatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

GeodeticLatLon(args...) = GeodeticLatLon{WGS84}(args...)

const LatLon = GeodeticLatLon

"""
    GeocentricLatLon(lat, lon)
    GeocentricLatLon{Datum}(lat, lon)

Geocentric latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
GeocentricLatLon(45, 45) # add default units
GeocentricLatLon(45u"°", 45u"°") # integers are converted converted to floats
GeocentricLatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
GeocentricLatLon(45.0u"°", 45.0u"°")
GeocentricLatLon{WGS84}(45.0u"°", 45.0u"°")
```
"""
struct GeocentricLatLon{Datum,D<:Deg} <: LatitudeLongitude{Datum}
  lat::D
  lon::D
  GeocentricLatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

GeocentricLatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = GeocentricLatLon{Datum}(promote(lat, lon)...)
GeocentricLatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = GeocentricLatLon{Datum}(rad2deg(lat), rad2deg(lon))
GeocentricLatLon{Datum}(lat::Number, lon::Number) where {Datum} =
  GeocentricLatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

GeocentricLatLon(args...) = GeocentricLatLon{WGS84}(args...)

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{GeocentricLatLon{Datum}}, (; lat, lon)::LatLon{Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan((1 - e²) * tan(lat)))
  GeocentricLatLon{Datum}(lat′ * u"°", lon)
end

function Base.convert(::Type{LatLon{Datum}}, (; lat, lon)::GeocentricLatLon{Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan(1 / (1 - e²) * tan(lat)))
  LatLon{Datum}(lat′ * u"°", lon)
end
