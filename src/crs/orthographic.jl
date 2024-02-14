# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Orthographic{lat₀,lon₀,S}

Orthographic CRS with latitude origin `lat₀` and longitude origin `lon₀` in degrees
and spherical mode `S` enabled or not.
"""
const Orthographic{lat₀,lon₀,S,M<:Met} = CRS{:Orthographic,@NamedTuple{x::M, y::M},WGS84,Tuple{lat₀,lon₀,S}}

Orthographic{lat₀,lon₀,S}(x::M, y::M) where {lat₀,lon₀,S,M<:Met} = Orthographic{lat₀,lon₀,S,float(M)}(x, y)
Orthographic{lat₀,lon₀,S}(x::Met, y::Met) where {lat₀,lon₀,S} = Orthographic{lat₀,lon₀,S}(promote(x, y)...)
Orthographic{lat₀,lon₀,S}(x::Len, y::Len) where {lat₀,lon₀,S} =
  Orthographic{lat₀,lon₀,S}(uconvert(u"m", x), uconvert(u"m", y))
Orthographic{lat₀,lon₀,S}(x::Number, y::Number) where {lat₀,lon₀,S} =
  Orthographic{lat₀,lon₀,S}(addunit(x, u"m"), addunit(y, u"m"))

"""
    OrthoNorth(x, y)

Orthographic North Pole coordinates in length units (default to meter).

## Examples

```julia
OrthoNorth(1, 1) # add default units
OrthoNorth(1u"m", 1u"m") # integers are converted converted to floats
OrthoNorth(1.0u"km", 1.0u"km") # length quantities are converted to meters
OrthoNorth(1.0u"m", 1.0u"m")
```
"""
const OrthoNorth = Orthographic{90.0u"°",0.0u"°",false}

"""
    OrthoSouth(x, y)

Orthographic South Pole coordinates in length units (default to meter).

## Examples

```julia
OrthoSouth(1, 1) # add default units
OrthoSouth(1u"m", 1u"m") # integers are converted converted to floats
OrthoSouth(1.0u"km", 1.0u"km") # length quantities are converted to meters
OrthoSouth(1.0u"m", 1.0u"m")
```
"""
const OrthoSouth = Orthographic{-90.0u"°",0.0u"°",false}

typealias(::Type{ESRI{102035}}) = Orthographic{90.0u"°",0.0u"°",true}

typealias(::Type{ESRI{102037}}) = Orthographic{-90.0u"°",0.0u"°",true}

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{Orthographic{lat₀,lon₀,false}}, coords::LatLon) where {lat₀,lon₀}
  🌎 = ellipsoid(Orthographic{lat₀,lon₀})
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(lon₀))
  ϕ₀ = oftype(ϕ, deg2rad(lat₀))
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  e² = oftype(l, eccentricity²(🌎))

  sinϕ = sin(ϕ)
  cosϕ = cos(ϕ)
  sinϕ₀ = sin(ϕ₀)
  cosϕ₀ = cos(ϕ₀)
  ν = 1 / sqrt(1 - e² * sinϕ^2)
  ν₀ = 1 / sqrt(1 - e² * sinϕ₀^2)

  x = a * (ν * cosϕ * sin(λ - λ₀))
  y = a * (ν * (sinϕ * cosϕ₀ - cosϕ * sinϕ₀ * cos(λ - λ₀)) + e² * (ν₀ * sinϕ₀ - ν * sinϕ) * cosϕ₀)

  Orthographic{lat₀,lon₀,false}(x * u"m", y * u"m")
end

function Base.convert(::Type{Orthographic{lat₀,lon₀,true}}, coords::LatLon) where {lat₀,lon₀}
  🌎 = ellipsoid(Orthographic{lat₀,lon₀})
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  λ₀ = oftype(λ, deg2rad(lon₀))
  ϕ₀ = oftype(ϕ, deg2rad(lat₀))
  l = ustrip(λ)
  a = oftype(l, ustrip(majoraxis(🌎)))

  cosϕ = cos(ϕ)
  x = a * cosϕ * sin(λ - λ₀)
  y = a * (sin(ϕ) * cos(ϕ₀) - cosϕ * sin(ϕ₀) * cos(λ - λ₀))

  Orthographic{lat₀,lon₀,true}(x * u"m", y * u"m")
end
