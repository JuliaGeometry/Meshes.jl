# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotation(rot)

Rotate geometry or mesh with rotation `rot` from ReferenceFrameRotations.jl. 

## Examples

```julia
Rotation(ReferenceFrameRotations.EulerAngleAxis(pi/4, [1, 0, 0]))
```

"""
struct Rotation{R} <: GeometricTransform
  rot::R
end

"""
    Rotation(axis, α; degrees=true)

Perform the rotation of angle `α` around the axis directed by the 
vector `axis`. The angle must be given in degrees if `degrees=true`, 
otherwise in radians.

"""
function Rotation(axis::Vector{Real}, α::Real; degrees::Bool = true)
  axisNorm = norm2(axis)
  if axisNorm == 0
    error("The `axis` vector is null.")
  end
  v = axis / axisNorm
  if degrees
    α = α * pi / 180
  end
  Rotation(EulerAngleAxis(α, v))
end

isrevertible(::Type{<:Rotation}) = true

function preprocess(transform::Rotation, object)
  rot = transform.rot
  local M, M⁻¹
  if typeof(rot) <: DCM
    M = rot
    M⁻¹ = inv_rotation(M)
  else
    M = convert(DCM, rot)
    M⁻¹ = convert(DCM, inv(rot))
  end
  M, M⁻¹
end

function applypoint(::Rotation, points, prep)
  M, _ = prep
  newpoints = [Point(M * p) for p in points]
  newpoints, prep
end

function revertpoint(::Rotation, newpoints, cache)
  _, M⁻¹ = cache
  [Point(M⁻¹ * p) for p in newpoints]
end

function reapplypoint(::Rotation, points, cache)
  M, _ = cache
  [Point(M * p) for p in points]
end
