# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

(::EuclideanDistance)(p₁::Point, p₂::Point) = norm(p₂ - p₁)

(d::EuclideanDistance)(p::Point, l::Line) = d(p, _nearestpoint(p, l))

(d::EuclideanDistance)(p::Point, s::Segment) = d(p, _nearestpoint(p, s))

function (d::EuclideanDistance)(l₁::Line, l₂::Line)
  λ₁, λ₂, r, rₐ = intersectparameters(l₁(0), l₁(1), l₂(0), l₂(1))
  if (r == rₐ == 2) || (r == rₐ == 1)  # lines intersect or are colinear
    zero(promote_type(lentype(l₁), lentype(l₂)))
  elseif (r == 1) && (rₐ == 2)  # lines are parallel
    d(l₁(0), l₂)
  else  # get distance between closest points on each line
    d(l₁(λ₁), l₂(λ₂))
  end
end

function (d::EuclideanDistance)(s₁::Segment, s₂::Segment)
  λ₁, λ₂, r, rₐ = intersectparameters(s₁(0), s₁(1), s₂(0), s₂(1))
  λ₁′ = clamp(λ₁, zero(λ₁), one(λ₁))
  λ₂′ = clamp(λ₂, zero(λ₂), one(λ₂))
  if (r == rₐ == 1) || (r == 1 && rₐ == 2)  # segments are colinear or parallel
    min(d(s₁(0), s₂), d(s₁(1), s₂))
  else  # get distance between closest points on each segment
    d(s₁(λ₁′), s₂(λ₂′))
  end
end

# exploit symmetry to get more methods for free
(d::EuclideanDistance)(g::Geometry, p::Point) = d(p, g)

# -----------------
# HELPER FUNCTIONS
# -----------------

function _nearestpoint(p::Point, l::Line)
  λ, a, b = _nearestparams(p, l)
  a + λ * (b - a)
end

function _nearestpoint(p::Point, s::Segment)
  λ, a, b = _nearestparams(p, s)
  λ′ = clamp(λ, zero(λ), one(λ))
  a + λ′ * (b - a)
end

function _nearestparams(p, l)
  a, b = l(0), l(1)
  u = p - a
  v = b - a
  λ = (u ⋅ v) / (v ⋅ v)
  λ, a, b
end
