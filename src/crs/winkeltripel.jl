# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Winkel{Datum,lat₁,M<:Met} <: CRS{Datum}
  x::M
  y::M
  Winkel{Datum,lat₁}(x::M, y::M) where {Datum,lat₁,M<:Met} = new{Datum,lat₁,float(M)}(x, y)
end

Winkel{Datum,lat₁}(x::Met, y::Met) where {Datum,lat₁} = Winkel{Datum,lat₁}(promote(x, y)...)
Winkel{Datum,lat₁}(x::Len, y::Len) where {Datum,lat₁} = Winkel{Datum,lat₁}(uconvert(u"m", x), uconvert(u"m", y))
Winkel{Datum,lat₁}(x::Number, y::Number) where {Datum,lat₁} = Winkel{Datum,lat₁}(addunit(x, u"m"), addunit(y, u"m"))

"""
    WinkelTripel{Datum}(x, y)

Winkel Tripel coordinates in length units (default to meter) with a given `Datum`.

## Examples

```julia
WinkelTripel{WGS84}(1, 1) # add default units
WinkelTripel{WGS84}(1u"m", 1u"m") # integers are converted converted to floats
WinkelTripel{WGS84}(1.0u"km", 1.0u"km") # length quantities are converted to meters
WinkelTripel{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54042](https://epsg.io/54042).
"""
const WinkelTripel{Datum} = Winkel{Datum,50.467u"°"}

typealias(::Type{ESRI{54042}}) = WinkelTripel{WGS84}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Winkel{Datum,lat₁}}, coords::LatLon{Datum}) where {Datum,lat₁}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕ₁ = oftype(ϕ, deg2rad(lat₁))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  α = acos(cos(ϕ) * cos(λ / 2))
  sincα = sinc(α / π) # unnormalized sinc

  x = a / 2 * (l * cos(ϕ₁) + (2cos(ϕ) * sin(λ / 2)) / sincα)
  y = a / 2 * (o + sin(ϕ) / sincα)

  Winkel{Datum,lat₁}(x * u"m", y * u"m")
end
