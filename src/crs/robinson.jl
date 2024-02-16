# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Robinson(x, y)
    Robinson{Datum}(x, y)

Robinson coordinates in length units (default to meter)
with a given `Datum` (default to `WGS84`).

## Examples

```julia
Robinson(1, 1) # add default units
Robinson(1u"m", 1u"m") # integers are converted converted to floats
Robinson(1.0u"km", 1.0u"km") # length quantities are converted to meters
Robinson(1.0u"m", 1.0u"m")
Robinson{WGS84}(1.0u"m", 1.0u"m")
```

See [ESRI:54030](https://epsg.io/54030).
"""
struct Robinson{Datum,M<:Met} <: CRS{Datum}
  x::M
  y::M
  Robinson{Datum}(x::M, y::M) where {Datum,M<:Met} = new{Datum,float(M)}(x, y)
end

typealias(::Type{ESRI{54030}}) = Robinson{WGS84}

Robinson(args...) = Robinson{WGS84}(args...)

Robinson{Datum}(x::Met, y::Met) where {Datum} = Robinson{Datum}(promote(x, y)...)
Robinson{Datum}(x::Len, y::Len) where {Datum} = Robinson{Datum}(uconvert(u"m", x), uconvert(u"m", y))
Robinson{Datum}(x::Number, y::Number) where {Datum} = Robinson{Datum}(addunit(x, u"m"), addunit(y, u"m"))

# ------------
# CONVERSIONS
# ------------

# Adapted from PROJ coordinate transformation software
# Initial PROJ 4.3 public domain code was put as Frank Warmerdam as copyright
# holder, but he didn't mean to imply he did the work. Essentially all work was
# done by Gerald Evenden.
# reference code: https://github.com/OSGeo/PROJ/blob/master/src/projections/robin.cpp

const _COEFSX = [
  (c₀=1.0f0, c₁=2.2199f-17, c₂=-7.15515f-5, c₃=3.1103f-6),
  (c₀=0.9986f0, c₁=-0.000482243f0, c₂=-2.4897f-5, c₃=-1.3309f-6),
  (c₀=0.9954f0, c₁=-0.00083103f0, c₂=-4.48605f-5, c₃=-9.86701f-7),
  (c₀=0.99f0, c₁=-0.00135364f0, c₂=-5.9661f-5, c₃=3.6777f-6),
  (c₀=0.9822f0, c₁=-0.00167442f0, c₂=-4.49547f-6, c₃=-5.72411f-6),
  (c₀=0.973f0, c₁=-0.00214868f0, c₂=-9.03571f-5, c₃=1.8736f-8),
  (c₀=0.96f0, c₁=-0.00305085f0, c₂=-9.00761f-5, c₃=1.64917f-6),
  (c₀=0.9427f0, c₁=-0.00382792f0, c₂=-6.53386f-5, c₃=-2.6154f-6),
  (c₀=0.9216f0, c₁=-0.00467746f0, c₂=-0.00010457f0, c₃=4.81243f-6),
  (c₀=0.8962f0, c₁=-0.00536223f0, c₂=-3.23831f-5, c₃=-5.43432f-6),
  (c₀=0.8679f0, c₁=-0.00609363f0, c₂=-0.000113898f0, c₃=3.32484f-6),
  (c₀=0.835f0, c₁=-0.00698325f0, c₂=-6.40253f-5, c₃=9.34959f-7),
  (c₀=0.7986f0, c₁=-0.00755338f0, c₂=-5.00009f-5, c₃=9.35324f-7),
  (c₀=0.7597f0, c₁=-0.00798324f0, c₂=-3.5971f-5, c₃=-2.27626f-6),
  (c₀=0.7186f0, c₁=-0.00851367f0, c₂=-7.01149f-5, c₃=-8.6303f-6),
  (c₀=0.6732f0, c₁=-0.00986209f0, c₂=-0.000199569f0, c₃=1.91974f-5),
  (c₀=0.6213f0, c₁=-0.010418f0, c₂=8.83923f-5, c₃=6.24051f-6),
  (c₀=0.5722f0, c₁=-0.00906601f0, c₂=0.000182f0, c₃=6.24051f-6),
  (c₀=0.5322f0, c₁=-0.00677797f0, c₂=0.000275608f0, c₃=6.24051f-6)
]

const _COEFSY = [
  (c₀=-5.20417f-18, c₁=0.0124f0, c₂=1.21431f-18, c₃=-8.45284f-11),
  (c₀=0.062f0, c₁=0.0124f0, c₂=-1.26793f-9, c₃=4.22642f-10),
  (c₀=0.124f0, c₁=0.0124f0, c₂=5.07171f-9, c₃=-1.60604f-9),
  (c₀=0.186f0, c₁=0.0123999f0, c₂=-1.90189f-8, c₃=6.00152f-9),
  (c₀=0.248f0, c₁=0.0124002f0, c₂=7.10039f-8, c₃=-2.24f-8),
  (c₀=0.31f0, c₁=0.0123992f0, c₂=-2.64997f-7, c₃=8.35986f-8),
  (c₀=0.372f0, c₁=0.0124029f0, c₂=9.88983f-7, c₃=-3.11994f-7),
  (c₀=0.434f0, c₁=0.0123893f0, c₂=-3.69093f-6, c₃=-4.35621f-7),
  (c₀=0.4958f0, c₁=0.0123198f0, c₂=-1.02252f-5, c₃=-3.45523f-7),
  (c₀=0.5571f0, c₁=0.0121916f0, c₂=-1.54081f-5, c₃=-5.82288f-7),
  (c₀=0.6176f0, c₁=0.0119938f0, c₂=-2.41424f-5, c₃=-5.25327f-7),
  (c₀=0.6769f0, c₁=0.011713f0, c₂=-3.20223f-5, c₃=-5.16405f-7),
  (c₀=0.7346f0, c₁=0.0113541f0, c₂=-3.97684f-5, c₃=-6.09052f-7),
  (c₀=0.7903f0, c₁=0.0109107f0, c₂=-4.89042f-5, c₃=-1.04739f-6),
  (c₀=0.8435f0, c₁=0.0103431f0, c₂=-6.4615f-5, c₃=-1.40374f-9),
  (c₀=0.8936f0, c₁=0.00969686f0, c₂=-6.4636f-5, c₃=-8.547f-6),
  (c₀=0.9394f0, c₁=0.00840947f0, c₂=-0.000192841f0, c₃=-4.2106f-6),
  (c₀=0.9761f0, c₁=0.00616527f0, c₂=-0.000256f0, c₃=-4.2106f-6),
  (c₀=1.0f0, c₁=0.00328947f0, c₂=-0.000319159f0, c₃=-4.2106f-6)
]

const _FXC = 0.8487
const _FYC = 1.3523
const _C₁ = 11.45915590261646417544
const _RC₁ = 0.08726646259971647884

function _V(C, z)
  c₀ = oftype(z, C.c₀)
  c₁ = oftype(z, C.c₁)
  c₂ = oftype(z, C.c₂)
  c₃ = oftype(z, C.c₃)
  c₀ + z * (c₁ + z * (c₂ + z * c₃))
end

function Base.convert(::Type{Robinson{Datum}}, coords::LatLon{Datum}) where {Datum}
  🌎 = ellipsoid(Datum)
  λ = deg2rad(coords.lon)
  ϕ = deg2rad(coords.lat)
  l = ustrip(λ)
  o = ustrip(ϕ)
  a = oftype(l, ustrip(majoraxis(🌎)))
  absϕ = abs(o)

  FXC = oftype(l, _FXC)
  FYC = oftype(l, _FYC)
  C₁ = oftype(l, _C₁)
  RC₁ = oftype(l, _RC₁)
  i = min(floor(Int, absϕ * C₁ + 1e-15), 18)
  z = rad2deg(absϕ - RC₁ * i)

  x = a * _V(_COEFSX[i + 1], z) * FXC * l
  y = a * _V(_COEFSY[i + 1], z) * FYC * sign(o)

  Robinson{Datum}(x * u"m", y * u"m")
end
