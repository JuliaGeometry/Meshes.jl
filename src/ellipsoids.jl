# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type RevolutionEllipsoid end

abstract type WGS84🌎 <: RevolutionEllipsoid end

const _WGS84 = let
  a = 6378137.0 * u"m"
  f⁻¹ = 298.257223563
  f = inv(f⁻¹)
  b = a * (1 - f)
  e² = (2 - f) / f⁻¹
  e = √e²
  (; a, b, e, e², f, f⁻¹)
end

majoraxis(::Type{WGS84🌎}) = _WGS84.a
minoraxis(::Type{WGS84🌎}) = _WGS84.b
eccentricity(::Type{WGS84🌎}) = _WGS84.e
eccentricity²(::Type{WGS84🌎}) = _WGS84.e
flattening(::Type{WGS84🌎}) = _WGS84.f
flattening⁻¹(::Type{WGS84🌎}) = _WGS84.f⁻¹
