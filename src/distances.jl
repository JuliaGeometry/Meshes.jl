# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# flip arguments so that points always come first
evaluate(d::PreMetric, g::Geometry, p::Point) = evaluate(d, p, g)

"""
    evaluate(distance::Euclidean, point, line)

Evaluate the Euclidean `distance` between `point` and `line`.
"""
function evaluate(::Euclidean, p::Point, l::Line)
  a, b = l(0), l(1)
  u = p - a
  v = b - a
  α = (u ⋅ v) / (v ⋅ v)
  norm(u - α * v)
end

"""
    evaluate(distance::Euclidean, line₁, line₂)

Evaluate the minimum Euclidean `distance` between `line₁` and `line₂`.
"""
function evaluate(d::Euclidean, l₁::Line, l₂::Line)
  λ₁, λ₂, r, rₐ = intersectparameters(l₁(0), l₁(1), l₂(0), l₂(1))

  if (r == rₐ == 2) || (r == rₐ == 1)  # lines intersect or are colinear
    return zero(result_type(d, lentype(l₁), lentype(l₂)))
  elseif (r == 1) && (rₐ == 2)  # lines are parallel
    return evaluate(d, l₁(0), l₂)
  else  # get distance between closest points on each line
    return evaluate(d, l₁(λ₁), l₂(λ₂))
  end
end

"""
    evaluate(distance::PreMetric, point₁, point₂)

Evaluate pre-metric `distance` between coordinates of `point₁` and `point₂`.
"""
function evaluate(d::PreMetric, p₁::Point, p₂::Point)
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  # TODO: maybe use cartvalues ​​with SVector
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
