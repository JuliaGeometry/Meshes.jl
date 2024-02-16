# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LatLon(lat, lon)
    LatLon{Datum}(lat, lon)

Latitude `lat ∈ [-90°,90°]` and longitude `lon ∈ [-180°,180°]` in angular units (default to degree)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
LatLon(45, 45) # add default units
LatLon(45u"°", 45u"°") # integers are converted converted to floats
LatLon((π/4)u"rad", (π/4)u"rad") # radians are converted to degrees
LatLon(45.0u"°", 45.0u"°")
LatLon{WGS84}(45.0u"°", 45.0u"°")
```

See [EPSG:4326](https://epsg.io/4326).
"""
struct LatLon{Datum,D<:Deg} <: CRS{Datum}
  lat::D
  lon::D
  LatLon{Datum}(lat::D, lon::D) where {Datum,D<:Deg} = new{Datum,float(D)}(lat, lon)
end

LatLon{Datum}(lat::Deg, lon::Deg) where {Datum} = LatLon{Datum}(promote(lat, lon)...)
LatLon{Datum}(lat::Rad, lon::Rad) where {Datum} = LatLon{Datum}(rad2deg(lat), rad2deg(lon))
LatLon{Datum}(lat::Number, lon::Number) where {Datum} = LatLon{Datum}(addunit(lat, u"°"), addunit(lon, u"°"))

LatLon(args...) = LatLon{WGS84}(args...)
