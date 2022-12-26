# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  RotateCoords(rot)

Rotate geometry or mesh with rotation `rot` from ReferenceFrameRotations.jl. 

## Examples

```julia
RotateCoords(EulerAngleAxis(pi/4, [1, 0, 0]))
```
"""
struct RotateCoords{R} <: StatelessGeometricTransform
  rot::R
end

"""
  RotateCoords(u, v)

Rotation mapping the axis directed by `u` to the axis directed by `v`. 
More precisely, it maps the plane passing through the origin with normal 
vector `u` to the plane passing through the origin with normal vector `v`.

## Examples

```julia
RotateCoords((1, 0, 0), (1, 1, 1))
```
"""
function RotateCoords(u::Tuple, v::Tuple)
  u⃗ = normalize(Vec(u))
  v⃗ = normalize(Vec(v))
  realpart = √((1 + u⃗ ⋅ v⃗)/2)
  imagpart = (u⃗ × v⃗) / 2realpart
  q = Quaternion(realpart, imagpart...)
  RotateCoords(q)
end

isrevertible(::Type{<:RotateCoords}) = true

function preprocess(transform::RotateCoords, object)
  rot = transform.rot
  convert.(DCM, (inv(rot), rot))
end

function applypoint(::RotateCoords, points, prep)
  R, _ = prep
  newpoints = [Point(R * coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::RotateCoords, newpoints, cache)
  _, R⁻¹ = cache
  [Point(R⁻¹ * coordinates(p)) for p in newpoints]
end