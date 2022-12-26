# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  Rotate(rot)

Rotate geometry or mesh with rotation `rot`
from ReferenceFrameRotations.jl. 

## Examples

```julia
Rotate(EulerAngleAxis(pi/4, [1, 0, 0]))
```
"""
struct Rotate{R} <: StatelessGeometricTransform
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
function Rotate(u::Vec, v::Vec)
  u⃗ = normalize(u)
  v⃗ = normalize(v)
  re = √((1 + u⃗ ⋅ v⃗)/2)
  im = (u⃗ × v⃗) / 2re
  Rotate(Quaternion(re, im...))
end

isrevertible(::Type{<:Rotate}) = true

function preprocess(transform::Rotate, object)
  rot = transform.rot
  convert.(DCM, (inv(rot), rot))
end

function applypoint(::Rotate, points, prep)
  R, _ = prep
  newpoints = [Point(R * coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::Rotate, newpoints, cache)
  _, R⁻¹ = cache
  [Point(R⁻¹ * coordinates(p)) for p in newpoints]
end