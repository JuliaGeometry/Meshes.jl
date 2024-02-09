# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mercator(x, y)

Mercator coordinates in length units (default to meter).

## Examples

```julia
Mercator(1, 1) # add default units
Mercator(1u"m", 1u"m") # integers are converted converted to floats
Mercator(1.0u"km", 1.0u"km") # length quantities are converted to meters
Mercator(1.0u"m", 1.0u"m")
```

See [EPSG:3395](https://epsg.io/3395).
"""
const Mercator{M<:Met} = CRS{EPSG{3395},@NamedTuple{x::M, y::M},WGS84}

typealias(::Type{EPSG{3395}}) = Mercator

Mercator(x::M, y::M) where {M<:Met} = Mercator{float(M)}(x, y)
Mercator(x::Met, y::Met) = Mercator(promote(x, y)...)
Mercator(x::Len, y::Len) = Mercator(uconvert(u"m", x), uconvert(u"m", y))
Mercator(x::Number, y::Number) = Mercator(addunit(x, u"m"), addunit(y, u"m"))

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Mercator}, coords::LatLon)
  🌎 = ellipsoid(coords)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  e = oftype(l, eccentricity(🌎))
  x = a * l
  y = a * (asinh(tan(ϕ)) - e * atanh(e * sin(ϕ)))
  Mercator(x * u"m", y * u"m")
end
