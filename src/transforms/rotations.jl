# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rotation(rotation)

Perform the rotation represented by `rotation`, a
`ReferenceFrameRotation` object from the `ReferenceFrameRotations` 
library.

## Examples

```julia
using ReferenceFrameRotations
Rotation(ReferenceFrameRotations.EulerAngleAxis(pi/4, [1, 0, 0]))
```

"""
struct Rotation{R} <: Meshes.GeometricTransform
  rotation::R
end

"""
    Rotation(axis, α; degrees=true)

Perform the rotation of angle `α` around the axis directed by the 
vector `axis`. The angle must be given in degrees if `degrees=true`, 
otherwise in radians.

"""
function Rotation(axis::Vector{Real}, α::Real; degrees::Bool = true)
  axisNorm = LinearAlgebra.norm2(axis)
  if axisNorm == 0
    error("The `axis` vector is null.")
  end
  v = axis / axisNorm
  if degrees
    α = α * pi / 180
  end
  R = ReferenceFrameRotations.EulerAngleAxis
  rotation = R(α, v)
  Rotation{R}(rotation)
end

isrevertible(::Type{<:Rotation}) = true

function preprocess(transform::Rotation, object)
  rotation = transform.rotation
  rtype = typeof(rotation)
  local M, invM
  if rtype <: ReferenceFrameRotations.DCM
    M = rotation
    invM = ReferenceFrameRotations.inv_rotation(M)
  elseif rtype <: ReferenceFrameRotations.EulerAngleAxis
    M = ReferenceFrameRotations.angleaxis_to_dcm(rotation)
    invM = ReferenceFrameRotations.inv(rotation)
  elseif rtype <: ReferenceFrameRotations.EulerAngles
    M = ReferenceFrameRotations.angle_to_dcm(rotation)
    invM = ReferenceFrameRotations.inv(rotation)
  elseif rtype <: ReferenceFrameRotations.Quaternion
    M = ReferenceFrameRotations.quat_to_dcm(rotation)
    invM = ReferenceFrameRotations.inv(rotation)
  end
  M, invM
end

function applypoint(::Rotation, points, prep)
  M, _ = prep
  newpoints = [Point(M * p) for p in points]
  newpoints, prep
end

function revertpoint(::Rotation, newpoints, cache)
  _, invM = cache
  [Point(invM * p) for p in newpoints]
end

function reapplypoint(::Rotation, points, cache)
  M, _ = cache
  [Point(M * p) for p in points]
end
