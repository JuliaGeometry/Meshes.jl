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