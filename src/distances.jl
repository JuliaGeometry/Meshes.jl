# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricDistance

A distance between points that is meaningful in physical settings.

Distances are functors, and are evaluated with a pair of points:

```julia
d = EuclideanDistance()

d(p₁, p₂)
```

The definition that is used depends on the manifold of the points,
which is available at compile time.
"""
abstract type GeometricDistance end

"""
    EuclideanDistance()

Length of the straight line connecting two points in the space
where they are embedded. For points on the ellipsoid (`🌐`) this
is the chord through the interior of the ellipsoid, which is
shorter than any path along the surface.

See also [`GeodesicDistance`](@ref).

## Examples

```julia
d = EuclideanDistance()

d(Point(0, 0), Point(1, 1))

d(Point(LatLon(0, 0)), Point(LatLon(0, 1)))
```
"""
struct EuclideanDistance <: GeometricDistance end

(::EuclideanDistance)(p₁::Point{M}, p₂::Point{M}) where {M} = norm(p₂ - p₁)

"""
    ManhattanDistance()

Sum of the absolute differences of the coordinates of two points
in Euclidean space (`𝔼`), also known as the taxicab distance.

## Examples

```julia
d = ManhattanDistance()

d(Point(0, 0), Point(1, 1))
```
"""
struct ManhattanDistance <: GeometricDistance end

(::ManhattanDistance)(p₁::Point{𝔼{Dim}}, p₂::Point{𝔼{Dim}}) where {Dim} = sum(abs, p₂ - p₁)

"""
    GeodesicDistance()

Length of the shortest path along the manifold of two points.
In Euclidean space (`𝔼`) this is the straight line connecting
them. On the ellipsoid (`🌐`) this is the geodesic of the
ellipsoid attached to the datum of the coordinate reference
system, and is computed with the series of Karney (2013).

See also [`EuclideanDistance`](@ref).

## Examples

```julia
d = GeodesicDistance()

d(Point(LatLon(0, 0)), Point(LatLon(0, 1)))
```

## References

* Karney, C. F. F. 2013. [Algorithms for geodesics]
  (https://doi.org/10.1007/s00190-012-0578-0)
"""
struct GeodesicDistance <: GeometricDistance end

(::GeodesicDistance)(p₁::Point{𝔼{Dim}}, p₂::Point{𝔼{Dim}}) where {Dim} = norm(p₂ - p₁)

function (::GeodesicDistance)(p₁::Point{🌐}, p₂::Point{🌐})
  q₁, q₂ = promote(p₁, p₂)
  c₁ = convert(manifoldcrs(q₁), coords(q₁))
  c₂ = convert(manifoldcrs(q₂), coords(q₂))

  # the manifold only tells us that the points lie on a sphere,
  # the ellipsoid itself comes from the datum of the coordinates
  🌎 = ellipsoid(datum(c₁))
  a = majoraxis(🌎)
  f = flattening(🌎)

  T = numtype(lentype(q₁))
  # the series of Karney need double precision to reach round-off
  S = promote_type(T, Float64)
  g = GeodesicEllipsoid(S(ustrip(a)), S(f))

  lat₁, lon₁ = S(ustrip(u"°", c₁.lat)), S(ustrip(u"°", c₁.lon))
  lat₂, lon₂ = S(ustrip(u"°", c₂.lat)), S(ustrip(u"°", c₂.lon))

  T(_geodesic(g, lat₁, lon₁, lat₂, lon₂)) * unit(a)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("distances/geodesic.jl")

# ---------------
# DISTANCES.JL
# ---------------

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
    evaluate(distance::Euclidean, point, segment)

Evaluate the Euclidean `distance` between `point` and `segment`.
"""
function evaluate(::Euclidean, p::Point, s::Segment)
  a, b = s(0), s(1)
  u = p - a
  v = b - a
  α = (u ⋅ v) / (v ⋅ v)
  α = clamp(α, zero(α), one(α))
  norm(u - α * v)
end

"""
    evaluate(distance::Euclidean, line₁, line₂)

Evaluate the minimum Euclidean `distance` between `line₁` and `line₂`.
"""
function evaluate(d::Euclidean, l₁::Line, l₂::Line)
  λ₁, λ₂, r, rₐ = intersectparameters(l₁(0), l₁(1), l₂(0), l₂(1))
  if (r == rₐ == 2) || (r == rₐ == 1)  # lines intersect or are colinear
    zero(result_type(d, lentype(l₁), lentype(l₂)))
  elseif (r == 1) && (rₐ == 2)  # lines are parallel
    evaluate(d, l₁(0), l₂)
  else  # get distance between closest points on each line
    evaluate(d, l₁(λ₁), l₂(λ₂))
  end
end

"""
    evaluate(distance::Euclidean, segment₁, segment₂)

Evaluate the minimum Euclidean `distance` between `segment₁` and `segment₂`.
"""
function evaluate(d::Euclidean, s₁::Segment, s₂::Segment)
  λ₁, λ₂, r, rₐ = intersectparameters(s₁(0), s₁(1), s₂(0), s₂(1))
  λ₁ = clamp(λ₁, zero(λ₁), one(λ₁))
  λ₂ = clamp(λ₂, zero(λ₂), one(λ₂))
  if (r == rₐ == 1) || (r == 1 && rₐ == 2)  # segments are colinear or parallel
    min(evaluate(d, s₁(0), s₂), evaluate(d, s₁(1), s₂))
  else  # get distance between closest points on each segment
    evaluate(d, s₁(λ₁), s₂(λ₂))
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
