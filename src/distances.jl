# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# flip arguments so that points always come first
evaluate(dist::PreMetric, g::Geometry, p::Point) = evaluate(dist, p, g)

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
function evaluate(dist::Euclidean, line₁::Line{Dim}, line₂::Line{Dim}) where {Dim}
  λ₁, λ₂, r, rₐ = intersectparameters(line₁(0), line₁(1), line₂(0), line₂(1))

  if (r == rₐ == 2) || (r == rₐ == 1)  # lines intersect or are colinear
    return zero(result_type(dist, lentype(line₁), lentype(line₂)))
  elseif (r == 1) && (rₐ == 2)  # lines are parallel
    return evaluate(dist, line₁(0), line₂)
  else  # get distance between closest points on each line
    return evaluate(dist, line₁(λ₁), line₂(λ₂))
  end
end

"""
    evaluate(distance::PreMetric, point₁, point₂)

Evaluate pre-metric `distance` between coordinates of `point₁` and `point₂`.
"""
function evaluate(dist::PreMetric, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim}
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  v₁ = ustrip.(u, to(p₁))
  v₂ = ustrip.(u, to(p₂))
  evaluate(dist, v₁, v₂) * u
end

# --------------
# SPECIAL CASES
# --------------

function evaluate(dist::Haversine, p₁::Point{Dim,<:LatLon}, p₂::Point{Dim,<:LatLon}) where {Dim}
  uᵣ = unit(dist.radius)
  u = uᵣ === NoUnits ? u"m" : uᵣ
  latlon₁ = coords(p₁)
  latlon₂ = coords(p₂)
  v₁ = SVector(latlon₁.lon, latlon₁.lat)
  v₂ = SVector(latlon₂.lon, latlon₂.lat)
  evaluate(dist, v₁, v₂) * u
end

evaluate(dist::Haversine, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim} =
  evaluate(dist, Point(convert(LatLon, coords(p₁))), Point(convert(LatLon, coords(p₂))))

function evaluate(dist::SphericalAngle, p₁::Point{Dim,<:LatLon}, p₂::Point{Dim,<:LatLon}) where {Dim}
  latlon₁ = coords(p₁)
  latlon₂ = coords(p₂)
  v₁ = SVector(deg2rad(latlon₁.lon), deg2rad(latlon₁.lat))
  v₂ = SVector(deg2rad(latlon₂.lon), deg2rad(latlon₂.lat))
  evaluate(dist, v₁, v₂) * u"rad"
end

evaluate(dist::SphericalAngle, p₁::Point{Dim}, p₂::Point{Dim}) where {Dim} =
  evaluate(dist, Point(convert(LatLon, coords(p₁))), Point(convert(LatLon, coords(p₂))))
