# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotate(rot)

Rotate geometry or mesh with rotation `rot`
from Rotations.jl.

## Examples

```julia
Rotate(one(RotXYZ{Float64}))  # Generate identity rotation
Rotate(AngleAxis(0.2, 1.0, 0.0, 0.0))  # Rotate 0.2 radians around X-axis
Rotate(rand(QuatRotation{Float64}))  # Generate random rotation
```
"""
struct Rotate{R<:Rotation} <: StatelessGeometricTransform
  rot::R
end

"""
    Rotate(u, v)

Rotation mapping the axis directed by `u` to the axis directed by `v`. 
More precisely, it maps the plane passing through the origin with normal 
vector `u` to the plane passing through the origin with normal vector `v`.

## Examples

```julia
Rotate(Vec(1, 0, 0), Vec(1, 1, 1))
```
"""
Rotate(u::Vec, v::Vec) = Rotate(rotation_between(u, v))

isrevertible(::Type{<:Rotate}) = true

Base.inv(transform::Rotate) = Rotate(inv(transform.rot))

function applypoint(transform::Rotate, points, prep)
  R = transform.rot
  newpoints = [Point(R * coordinates(p)) for p in points]
  newpoints, inv(R)
end

function revertpoint(::Rotate, newpoints, cache)
  R⁻¹ = cache
  [Point(R⁻¹ * coordinates(p)) for p in newpoints]
end

# --------------
# SPECIAL CASES
# --------------

function apply(transform::Rotate, plane::Plane)
  R = transform.rot
  o = plane(0, 0)
  n = normal(plane)
  Plane(o, R * n), inv(R)
end

function revert(::Rotate, plane::Plane, cache)
  R⁻¹ = cache
  o = plane(0, 0)
  n = normal(plane)
  Plane(o, R⁻¹ * n)
end

function apply(transform::Rotate, cylinder::Cylinder)
  bot, top = _rcylinder(transform, cylinder)
  Cylinder(bot, top, radius(cylinder)), nothing
end

function revert(transform::Rotate, cylinder::Cylinder, cache)
  bot, top = _rcylinder(inv(transform), cylinder)
  Cylinder(bot, top, radius(cylinder))
end

function _rcylinder(transform, cylinder)
  b = bottom(cylinder)
  t = top(cylinder)
  transform(b), transform(t)
end
