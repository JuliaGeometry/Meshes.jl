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