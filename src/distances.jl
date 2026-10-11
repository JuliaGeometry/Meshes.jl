# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricDistance

A distance between points in geometries.
"""
abstract type GeometricDistance end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("distances/euclidean.jl")
include("distances/geodesic.jl")

# -------------
# DISTANCES.JL
# -------------

"""
    evaluate(distance::PreMetric, point₁, point₂)

Evaluate pre-metric `distance` between coordinates of `point₁` and `point₂`.
"""
function evaluate(d::PreMetric, p₁::Point, p₂::Point)
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  v₁ = ustrip.(u, to(p₁))
  v₂ = ustrip.(u, to(p₂))
  evaluate(d, v₁, v₂) * u
end

# --------------
# SPECIAL CASES
# --------------

evaluate(d::Haversine, p₁::Point, p₂::Point) = _evaluate(d, coords(p₁), coords(p₂))

function _evaluate(d::Haversine, coords₁::LatLon, coords₂::LatLon)
  uᵣ = unit(d.radius)
  # add default unit if necessary
  u = uᵣ === NoUnits ? u"m" : NoUnits
  v₁ = SVector(coords₁.lon, coords₁.lat)
  v₂ = SVector(coords₂.lon, coords₂.lat)
  evaluate(d, v₁, v₂) * u
end

_evaluate(d::Haversine, coords₁::CRS, coords₂::CRS) = _evaluate(d, convert(LatLon, coords₁), convert(LatLon, coords₂))

evaluate(d::SphericalAngle, p₁::Point, p₂::Point) = _evaluate(d, coords(p₁), coords(p₂))

function _evaluate(d::SphericalAngle, coords₁::LatLon, coords₂::LatLon)
  v₁ = SVector(deg2rad(coords₁.lon), deg2rad(coords₁.lat))
  v₂ = SVector(deg2rad(coords₂.lon), deg2rad(coords₂.lat))
  evaluate(d, v₁, v₂) * u"rad"
end

_evaluate(d::SphericalAngle, coords₁::CRS, coords₂::CRS) =
  _evaluate(d, convert(LatLon, coords₁), convert(LatLon, coords₂))
