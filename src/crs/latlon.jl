# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
const LatLon{D<:Deg} = CRS{EPSG{4326},@NamedTuple{lat::D, lon::D},WGS84,NoParams}

typealias(::Type{EPSG{4326}}) = LatLon

LatLon(lat::D, lon::D) where {D<:Deg} = LatLon{float(D)}(lat, lon)
LatLon(lat::Deg, lon::Deg) = LatLon(promote(lat, lon)...)
LatLon(lat::Rad, lon::Rad) = LatLon(rad2deg(lat), rad2deg(lon))
LatLon(lat::Number, lon::Number) = LatLon(addunit(lat, u"°"), addunit(lon, u"°"))
