@enum LatLonKind begin
  Geodetic
  Geocentric
end

"""
    LatLon(lat, lon)
    LatLon{Kind}(lat, lon)
    LatLon{Kind,Datum}(lat, lon)

Latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a specified `Kind` (`Geodetic`, the default, or `Geocentric`) and a given `Datum` (default to `WGS84`).

## Examples

```julia
LatLon(45, 45) # add default units
LatLon(45u"°", 45u"°") # integers are converted converted to floats
LatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
LatLon(45.0u"°", 45.0u"°")
LatLon{Geocentric}(45.0u"°", 45.0u"°")
LatLon{Geodetic,WGS84}(45.0u"°", 45.0u"°")
```

See [EPSG:4326](https://epsg.io/4326).
"""
struct LatLon{Kind,Datum,D<:Deg} <: CRS{Datum}
  lat::D
  lon::D
  LatLon{Kind,Datum}(lat::D, lon::D) where {Kind,Datum,D<:Deg} = new{Kind,Datum,float(D)}(lat, lon)
end

LatLon{Kind,Datum}(lat::Deg, lon::Deg) where {Kind,Datum} = LatLon{Kind,Datum}(promote(lat, lon)...)
LatLon{Kind,Datum}(lat::Rad, lon::Rad) where {Kind,Datum} = LatLon{Kind,Datum}(rad2deg(lat), rad2deg(lon))
LatLon{Kind,Datum}(lat::Number, lon::Number) where {Kind,Datum} =
  LatLon{Kind,Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

LatLon{Kind}(args...) where {Kind} = LatLon{Kind,WGS84}(args...)

LatLon(args...) = LatLon{Geodetic}(args...)

Base.summary(io::IO, ::LatLon{Kind,Datum}) where {Kind,Datum} = print(io, "($Kind) LatLon{$Datum} coordinates")

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{LatLon{Geocentric,Datum}}, (; lat, lon)::LatLon{Geodetic,Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan((1 - e²) * tan(lat)))
  LatLon{Geocentric,Datum}(lat′ * u"°", lon)
end

function Base.convert(::Type{LatLon{Geodetic,Datum}}, (; lat, lon)::LatLon{Geocentric,Datum}) where {Datum}
  l = ustrip(lat)
  e² = oftype(l, eccentricity²(ellipsoid(Datum)))
  lat′ = rad2deg(atan(1 / (1 - e²) * tan(lat)))
  LatLon{Geodetic,Datum}(lat′ * u"°", lon)
end
