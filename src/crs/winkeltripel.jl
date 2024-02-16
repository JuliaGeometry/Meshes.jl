# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Winkel{lat₁,Datum,M<:Met} <: CRS{Datum}
  x::M
  y::M
  Winkel{lat₁,Datum}(x::M, y::M) where {lat₁,Datum,M<:Met} = new{lat₁,Datum,float(M)}(x, y)
end

Winkel{lat₁}(args...) where {lat₁} = Winkel{lat₁,WGS84}(args...)

Winkel{lat₁,Datum}(x::Met, y::Met) where {lat₁,Datum} = Winkel{lat₁,Datum}(promote(x, y)...)
Winkel{lat₁,Datum}(x::Len, y::Len) where {lat₁,Datum} = Winkel{lat₁,Datum}(uconvert(u"m", x), uconvert(u"m", y))
Winkel{lat₁,Datum}(x::Number, y::Number) where {lat₁,Datum} = Winkel{lat₁,Datum}(addunit(x, u"m"), addunit(y, u"m"))

"""
    WinkelTripel(x, y)
    WinkelTripel{Datum}(x, y)

Winkel Tripel coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
WinkelTripel(1, 1) # add default units
WinkelTripel(1u"m", 1u"m") # integers are converted converted to floats
WinkelTripel(1.0u"km", 1.0u"km") # length quantities are converted to meters
WinkelTripel(1.0u"m", 1.0u"m")
WinkelTripel{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54042](https://epsg.io/54042).
"""
const WinkelTripel{Datum} = Winkel{50.467u"°",Datum}

typealias(::Type{ESRI{54042}}) = WinkelTripel{WGS84}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Winkel{lat₁,Datum}}, coords::LatLon{Datum}) where {lat₁,Datum}
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

  Winkel{lat₁,Datum}(x * u"m", y * u"m")
end
