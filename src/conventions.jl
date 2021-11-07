# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RotationConvention

Convention used to specify rotations in 2D and 3D spaces.
"""
abstract type RotationConvention end

"""
    axesseq(convention) -> Symbol

Sequence of three axes by which the rotations are made, e.g. `:ZXZ`.
"""
function axesseq(::Type{<:RotationConvention}) end

"""
    orientation(convention) -> Tuple{Symbol,Symbol,Symbol}

Informs for each of the three rotations if it is clockwise (`:CW`) or
counter-clockwise (`:CCW`). Adopts right-hand rule: orientation defined
looking towards the negative direction of the axis.
"""
function orientation(::Type{<:RotationConvention}) end

"""
    angleunits(convention) -> Symbol

Informs the angle units of the rotations in radians (`:RAD`) or in degrees (`:DEG`).
"""
function angleunits(::Type{<:RotationConvention}) end

"""
    mainaxis(convention) -> Symbol

Informs if the main axis of the rotation is `:X` or `:Y`.
"""
function mainaxis(::Type{<:RotationConvention}) end

"""
    isextrinsic(convention) -> Bool

Informs if the rotation is extrinsic or intrinsic.
"""
function isextrinsic(::Type{<:RotationConvention}) end

"""
    rotmat(angles, convention)

Return the direction cosine matrix of rotation for the given `angles`
according to the rotation `convention`.
"""
function rotmat(angles, convention::Type{<:RotationConvention})
  nangles = length(angles)
  @assert nangles ∈ [1,3] "invalid number of angles"

  Dim = nangles == 1 ? 2 : 3

  # convert to radian if necessary
  θs = angleunits(convention) == :DEG ? deg2rad.(angles) : angles

  # invert sign if necessary
  _0 = zero(eltype(θs))
  Dim == 2 && (θs = [θs[1], _0, _0])
  intr = (orientation(convention) .== :CW) .& !isextrinsic(convention)
  extr = (orientation(convention) .== :CCW) .& isextrinsic(convention)
  inds = collect(intr .| extr)
  θs[inds] *= -1

  # rotation matrix
  R = angle_to_dcm(θs..., axesseq(convention))[SOneTo(Dim),SOneTo(Dim)]

  isextrinsic(convention) ? R' : R
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("conventions/TaitBryanExtr.jl")
include("conventions/TaitBryanIntr.jl")
include("conventions/EulerExtr.jl")
include("conventions/EulerIntr.jl")
include("conventions/GSLIB.jl")
include("conventions/Leapfrog.jl")
include("conventions/Datamine.jl")
