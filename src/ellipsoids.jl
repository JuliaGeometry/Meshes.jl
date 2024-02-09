# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RevolutionEllipsoid

Parent type of all revolution ellipsoids.
"""
abstract type RevolutionEllipsoid end

"""
    majoraxis(E)

Returns the semi-major axis of the ellipsoid type `E`.
"""
function majoraxis end

"""
    minoraxis(E)

Returns the semi-minor axis of the ellipsoid type `E`.
"""
function minoraxis end

"""
    eccentricity(E)

Returns the eccentricity of the ellipsoid type `E`.
"""
function eccentricity end

"""
    eccentricity²(E)

Returns the eccentricity squared of the ellipsoid type `E`.
"""
function eccentricity² end

"""
    flattening(E)

Returns the flattening of the ellipsoid type `E`.
"""
function flattening end

"""
    flattening⁻¹(E)

Returns the inverse flattening of the ellipsoid type `E`.
"""
function flattening⁻¹ end

abstract type WGS84🌎 <: RevolutionEllipsoid end

const _WGS84 = let
  a = 6378137.0u"m"
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
eccentricity²(::Type{WGS84🌎}) = _WGS84.e²
flattening(::Type{WGS84🌎}) = _WGS84.f
flattening⁻¹(::Type{WGS84🌎}) = _WGS84.f⁻¹

abstract type WIII🌎 <: RevolutionEllipsoid end

majoraxis(::Type{WIII🌎}) = 6371000.0u"m"
minoraxis(::Type{WIII🌎}) = 6371000.0u"m"
eccentricity(::Type{WIII🌎}) = 0.0
eccentricity²(::Type{WIII🌎}) = 0.0
flattening(::Type{WIII🌎}) = 0.0
flattening⁻¹(::Type{WIII🌎}) = Inf
