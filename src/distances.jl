# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# flip arguments so that points always come first
evaluate(d::PreMetric, g::Geometry, p::Point) = evaluate(d, p, g)

"""
    evaluate(distance::Euclidean, point, line)

Evaluate the Euclidean `distance` between `point` and `line`.
"""
function evaluate(::Euclidean, p::Point{Dim}, l::Line{Dim}) where {Dim}
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
function evaluate(d::Euclidean, l₁::Line{Dim}, l₂::Line{Dim}) where {Dim}
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
function evaluate(d::PreMetric, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim}
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

function evaluate(d::Haversine, p₁::Point{Dim,<:LatLon}, p₂::Point{Dim,<:LatLon}) where {Dim}
  uᵣ = unit(d.radius)
  # add default unit if necessary
  u = uᵣ === NoUnits ? u"m" : NoUnits
  latlon₁ = coords(p₁)
  latlon₂ = coords(p₂)
  v₁ = SVector(latlon₁.lon, latlon₁.lat)
  v₂ = SVector(latlon₂.lon, latlon₂.lat)
  evaluate(d, v₁, v₂) * u
end

evaluate(d::Haversine, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim} =
  evaluate(d, Point(convert(LatLon, coords(p₁))), Point(convert(LatLon, coords(p₂))))

function evaluate(d::SphericalAngle, p₁::Point{Dim,<:LatLon}, p₂::Point{Dim,<:LatLon}) where {Dim}
  latlon₁ = coords(p₁)
  latlon₂ = coords(p₂)
  v₁ = SVector(deg2rad(latlon₁.lon), deg2rad(latlon₁.lat))
  v₂ = SVector(deg2rad(latlon₂.lon), deg2rad(latlon₂.lat))
  evaluate(d, v₁, v₂) * u"rad"
end

evaluate(d::SphericalAngle, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim} =
  evaluate(d, Point(convert(LatLon, coords(p₁))), Point(convert(LatLon, coords(p₂))))
