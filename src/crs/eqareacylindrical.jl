# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EqualAreaCylindrical{latₜₛ}

Equal Area Cylindrical CRS with latitude of true scale `latₜₛ` in degrees.
"""
const EqualAreaCylindrical{latₜₛ,M<:Met} = CRS{:EAC,@NamedTuple{x::M, y::M},WGS84,latₜₛ}

EqualAreaCylindrical{latₜₛ}(x::M, y::M) where {latₜₛ,M<:Met} = EqualAreaCylindrical{latₜₛ,float(M)}(x, y)
EqualAreaCylindrical{latₜₛ}(x::Met, y::Met) where {latₜₛ} = EqualAreaCylindrical{latₜₛ}(promote(x, y)...)
EqualAreaCylindrical{latₜₛ}(x::Len, y::Len) where {latₜₛ} =
  EqualAreaCylindrical{latₜₛ}(uconvert(u"m", x), uconvert(u"m", y))
EqualAreaCylindrical{latₜₛ}(x::Number, y::Number) where {latₜₛ} =
  EqualAreaCylindrical{latₜₛ}(addunit(x, u"m"), addunit(y, u"m"))

"""
    Lambert(x, y)

Lambert cylindrical equal-area coordinates in length units (default to meter).

## Examples

```julia
Lambert(1, 1) # add default units
Lambert(1u"m", 1u"m") # integers are converted converted to floats
Lambert(1.0u"km", 1.0u"km") # length quantities are converted to meters
Lambert(1.0u"m", 1.0u"m")
```

See [ESRI:54034](https://epsg.io/54034).
"""
const Lambert = EqualAreaCylindrical{0.0u"°"}

typealias(::Type{ESRI{54034}}) = Lambert

"""
    Behrmann(x, y)

Behrmann coordinates in length units (default to meter).

## Examples

```julia
Lambert(1, 1) # add default units
Lambert(1u"m", 1u"m") # integers are converted converted to floats
Lambert(1.0u"km", 1.0u"km") # length quantities are converted to meters
Lambert(1.0u"m", 1.0u"m")
```

See [ESRI:54017](https://epsg.io/54017).
"""
const Behrmann = EqualAreaCylindrical{30.0u"°"}

typealias(::Type{ESRI{54017}}) = Behrmann

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{EqualAreaCylindrical{latₜₛ}}, coords::LatLon) where {latₜₛ}
  dat = datum(EqualAreaCylindrical{latₜₛ})
  🌎 = ellipsoid(dat)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(longitudeₒ(dat)))
  ϕₜₛ = oftype(ϕ, deg2rad(latₜₛ))
  l = ustrip(λ)
  l₀ = ustrip(λ₀)
  a = oftype(l, ustrip(majoraxis(🌎)))
  e = oftype(l, eccentricity(🌎))
  e² = oftype(l, eccentricity²(🌎))

  k₀ = cos(ϕₜₛ) / sqrt(1 - e² * sin(ϕₜₛ)^2)
  sinϕ = sin(ϕ)
  esinϕ = e * sinϕ
  q = (1 - e²) * (sinϕ / (1 - esinϕ^2) - (1 / 2e) * log((1 - esinϕ) / (1 + esinϕ)))

  x = a * k₀ * (l - l₀)
  y = a * q / 2k₀

  EqualAreaCylindrical{latₜₛ}(x * u"m", y * u"m")
end
