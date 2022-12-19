# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotation(rot)

Rotate geometry or mesh with rotation `rot` from ReferenceFrameRotations.jl. 

## Examples

```julia
Rotation(EulerAngleAxis(pi/4, [1, 0, 0]))
```
"""
struct Rotation{R} <: GeometricTransform
  rot::R
end

isrevertible(::Type{<:Rotation}) = true

function preprocess(transform::Rotation, object)
  rot = transform.rot
  R, R⁻¹ = rot, inv_rotation(rot)
  convert.(DCM, (R, R⁻¹))
end

function applypoint(::Rotation, points, prep)
  M, _ = prep
  newpoints = [Point(M * coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::Rotation, newpoints, cache)
  _, M⁻¹ = cache
  [Point(M⁻¹ * coordinates(p)) for p in newpoints]
end

function reapplypoint(::Rotation, points, cache)
  M, _ = cache
  [Point(M * coordinates(p)) for p in points]
end
