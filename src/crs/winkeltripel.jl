# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Winkel{ID,lat₁,M<:Met} = CRS{ID,@NamedTuple{x::M, y::M},WIII,lat₁}

Winkel{ID,lat₁}(x::M, y::M) where {ID,lat₁,M<:Met} = Winkel{ID,lat₁,float(M)}(x, y)
Winkel{ID,lat₁}(x::Met, y::Met) where {ID,lat₁} = Winkel{ID,lat₁}(promote(x, y)...)
Winkel{ID,lat₁}(x::Len, y::Len) where {ID,lat₁} = Winkel{ID,lat₁}(uconvert(u"m", x), uconvert(u"m", y))
Winkel{ID,lat₁}(x::Number, y::Number) where {ID,lat₁} = Winkel{ID,lat₁}(addunit(x, u"m"), addunit(y, u"m"))

"""
    WinkelTripel(x, y)

Winkel Tripel coordinates in length units (default to meter).

## Examples

```julia
WinkelTripel(1, 1) # add default units
WinkelTripel(1u"m", 1u"m") # integers are converted converted to floats
WinkelTripel(1.0u"km", 1.0u"km") # length quantities are converted to meters
WinkelTripel(1.0u"m", 1.0u"m")
```

See [ESRI:3395](https://epsg.io/53042).
"""
const WinkelTripel = Winkel{ESRI{53042},50.467u"°"}

typealias(::Type{ESRI{53042}}) = WinkelTripel

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Winkel{ID,lat₁}}, coords::LatLon) where {ID,lat₁}
  🌎 = ellipsoid(Winkel{ID,lat₁})
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  ϕ₁ = oftype(ϕ, deg2rad(lat₁))
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  α = acos(cos(ϕ) * cos(λ / 2))
  sincα = sin(α) / α # unnormalized sinc
  x = a / 2 * (l * cos(ϕ₁) + (2cos(ϕ) * sin(λ / 2)) / sincα)
  y = a / 2 * (o + sin(ϕ) / sincα)
  Winkel{ID,lat₁}(x * u"m", y * u"m")
end
