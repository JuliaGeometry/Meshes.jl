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
  # this method is a combination of a Gram-Schmidt solve
  # and the clamp and reproject method of Ericson,
  # "Real-Time Collision Detection" (2005), §5.1.9
  # http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  v₁ = s₁(1) - s₁(0)
  v₂ = s₂(1) - s₂(0)
  r = s₂(0) - s₁(0)
  a = v₁ ⋅ v₁
  b = v₁ ⋅ v₂
  c = v₂ ⋅ v₂
  f = v₁ ⋅ r
  e = v₂ ⋅ r
  a⁻¹ = inv(a)
  av₂ = b * a⁻¹
  ar = f * a⁻¹
  v₂p = fma.(-av₂, v₁, v₂) # component of v₂ ⟂ to v₁
  rp = fma.(-ar, v₁, r) # component of r ⟂ to v₁
  # in the plane ⟂ to v₁ find the nearest point on s₂
  λ₂ = _clamp01nan(-(rp ⋅ v₂p) / (v₂p ⋅ v₂p))
  # then find the nearest point on s₁ to that (we can reuse the divide by a)
  λ₁ = _clamp01nan(fma(λ₂, av₂, ar))
  # then find the nearest point on s₂ to that
  λ₂′ = _clamp01nan(fma(λ₁, b, -e) / c)
  norm(fma.(λ₂′, v₂, fma.(-λ₁, v₁, r)))
end

# exploit symmetry to get more methods for free
(d::EuclideanDistance)(g::Geometry, p::Point) = d(p, g)

# -----------------
# HELPER FUNCTIONS
# -----------------

# NaN λ becomes zero to handle degenerate segments
function _clamp01nan(λ)
  λ′ = ifelse(λ > zero(λ), λ, zero(λ))
  ifelse(λ′ < one(λ), λ′, one(λ))
end

function _nearestpoint(p::Point, l::Line)
  λ, a, b = _nearestparams(p, l)
  a + λ * (b - a)
end

function _nearestpoint(p::Point, s::Segment)
  λ, a, b = _nearestparams(p, s)
  λ′ = _clamp01nan(λ)
  a + λ′ * (b - a)
end

function _nearestparams(p, l)
  a, b = l(0), l(1)
  u = p - a
  v = b - a
  λ = (u ⋅ v) / (v ⋅ v)
  λ, a, b
end
